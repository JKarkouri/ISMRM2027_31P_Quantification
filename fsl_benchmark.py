#!/usr/bin/env python
"""FSL-MRS arm of the 31P NAD benchmark.

Fits the simulated spectra exported by fsl_export.m with FSL-MRS (Newton, Voigt
model), with the OLD and the Ren-corrected (NEW) basis, and writes the results
in the same layout as the LCModel arm so that bench_compare.m can merge them.

    python fsl_benchmark.py /mnt/c/.../fsl_bench/fsl_input.mat --workers 8
    python fsl_benchmark.py fsl_input.mat --configs "FSL-NEW poly2" --scen 4 5 6

Steps
  1. Calibration / convention check per configuration: a near-noiseless
     spectrum made from that basis (true amplitudes, +LW_CAL Hz) is fitted; the
     relative amplitudes must come back (rms log error < 0.05), otherwise the
     conjugated variants are tried and the script stops if none works. The
     absolute FSL-MRS scale is calibrated on the same fit (as for LCModel).
  2. Every spectrum x configuration is fitted (in parallel over spectra).
  3. EST_L [nS nQ nK nC nR nN], cfgL, T_FIT (s per fit) -> fsl_results.mat.

Use the Python that has fsl_mrs (e.g. $FSLDIR/bin/python, or your fsl_mrs
conda env). Check:  python -c "import fsl_mrs; print(fsl_mrs.__version__)"
"""
import os
for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')          # one thread per worker process
import argparse
import time
import warnings
from multiprocessing import Pool

import numpy as np
from scipy.io import loadmat, savemat

# name: (basis, baseline, metabolite grouping)
#   baseline: FSL-MRS string, 'poly, 2' is the fsl_mrs default; 'spline, moderate'
#   is the penalised spline (closer to the spline baselines of the other fitters,
#   but much slower).  grouping: None = one group (fsl_mrs default), 'grp' = by
#   resonance family (own shift / width per group).
CONFIGS = {
    'FSL-OLD poly2':  ('OLD', 'poly, 2', None),
    'FSL-NEW poly2':  ('NEW', 'poly, 2', None),
    'FSL-OLD spline': ('OLD', 'spline, moderate', None),
    'FSL-NEW spline': ('NEW', 'spline, moderate', None),
    # grouped: own shift / width per resonance family, as in the other two fitters
    'FSL-OLD poly2 grp':  ('OLD', 'poly, 2', 'grp'),
    'FSL-NEW poly2 grp':  ('NEW', 'poly, 2', 'grp'),
    'FSL-OLD spline grp': ('OLD', 'spline, moderate', 'grp'),
    'FSL-NEW spline grp': ('NEW', 'spline, moderate', 'grp'),
    # as grp, but NAD/UDP share the alpha-ATP group (FSL-MRS has no priors, so this is
    # its version of anchoring NAD to alpha-ATP)
    'FSL-OLD spline grpA': ('OLD', 'spline, moderate', 'grpA'),
    'FSL-NEW spline grpA': ('NEW', 'spline, moderate', 'grpA'),
}
DEFAULT_CONFIGS = ['FSL-OLD poly2', 'FSL-NEW poly2', 'FSL-OLD spline', 'FSL-NEW spline']
# same families as the in-house model: PME, PDE, Pi, NAD/UDP, and each ATP site separately
GROUPS = {'PCr': 0, 'Pi_intra': 1, 'Pi_extra': 1, 'PE': 2, 'PC': 2, 'GPC': 3, 'GPE': 3,
          'gATP': 4, 'aATP': 5, 'bATP': 6, 'NADH': 7, 'NADplus': 7, 'UDPG': 7,
          'UDPGal': 7, 'UDPGlc': 7, 'UDPGalNAc': 7, 'UDPGlcNAc': 7}
PPMLIM = (-20.0, 10.0)
CAL_METS = ['PCr', 'Pi_intra', 'gATP', 'aATP', 'GPC', 'PE', 'GPE']

WARM_START = True   # spline fits start from a poly-2 fit (--no-warm: FSL-MRS default init)
CACHE_LAMBDA = True  # reuse the spline penalty (identical result; --no-cache: FSL-MRS as shipped)
_LAM = {}


def _install_lambda_cache():
    """Memoise fsl_mrs' lambda_from_ed. The spline penalty depends only on the
    requested effective dimension and the baseline basis, which are the same for
    every spectrum here, but FSL-MRS recomputes it (twice) in every fit."""
    import hashlib
    import fsl_mrs.utils.baseline as bl
    if getattr(bl.lambda_from_ed, '_cached', False):
        return
    orig = bl.lambda_from_ed

    def cached(target_ed, basis):
        key = (float(target_ed), basis.shape, hashlib.sha1(np.ascontiguousarray(basis).tobytes()).hexdigest())
        if key not in _LAM:
            _LAM[key] = orig(target_ed, basis)
        return _LAM[key]
    cached._cached = True
    bl.lambda_from_ed = cached


D = {}   # data shared with the worker processes (fork)


def cellstr(x):
    return [str(np.squeeze(v)) for v in np.ravel(x)]


def load_input(path):
    m = loadmat(path)
    dims = [int(v) for v in np.ravel(m['dims'])]
    D['FID'] = np.reshape(m['FID'], dims, order='F') if m['FID'].ndim < 5 else m['FID']
    D['dt'] = float(m['dt'].squeeze())
    D['B0'] = float(m['B0'].squeeze())
    D['LW_CAL'] = float(m['LW_CAL'].squeeze())
    D['REP'] = cellstr(m['REP'])
    for b in ('OLD', 'NEW'):
        D['F_' + b] = np.asarray(m['F_' + b], complex)
        D['mets_' + b] = cellstr(m['mets_' + b])
        D['W_' + b] = np.asarray(m['W_' + b], float)
        D['aCal_' + b] = np.ravel(m['aCal_' + b]).astype(float)


def groups_for(mets, mode):
    if mode is None:
        return None
    tab = dict(GROUPS)
    if mode == 'grpA':
        tab.update({n: GROUPS['aATP'] for n, v in GROUPS.items() if v == GROUPS['NADH']})
    g = [tab.get(n, max(tab.values()) + 1) for n in mets]
    used = sorted(set(g))                          # renumber 0..G-1
    return [used.index(v) for v in g]


def fit_one(fid, bname, baseline, gmode, conj_fid=False, conj_basis=False, plot=None, return_res=False):
    """Raw FSL-MRS concentrations divided by the FID scaling (linear in the data)."""
    from fsl_mrs.core.mrs import MRS
    from fsl_mrs.utils.fitting import fit_FSLModel
    if CACHE_LAMBDA:
        _install_lambda_cache()
    F, mets = D['F_' + bname], D['mets_' + bname]
    hdr = [{'dwelltime': D['dt'], 'bandwidth': 1 / D['dt'], 'centralFrequency': D['B0'],
            'fwhm': 5.0}] * len(mets)
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        mrs = MRS(FID=fid.copy(), cf=D['B0'], bw=1 / D['dt'], nucleus='31P',
                  basis=F.copy(), names=mets, basis_hdr=hdr)
        if conj_fid:
            mrs.conj_FID = True
        if conj_basis:
            mrs.conj_Basis = True
        mrs.rescaleForFitting()
        grp = groups_for(mets, gmode)
        x0 = None
        g = 1 if grp is None else max(grp) + 1
        if WARM_START and (baseline.startswith('spline') or g > 1):
            # Stage 1: robust single-group polynomial fit. Its amplitudes, phase and
            # lineshape (copied to every group) start the final fit, whose own
            # baseline starts at zero. Without this, the spline and grouped fits
            # can get stuck (baseline absorbing peaks, metabolites set to zero).
            from fsl_mrs.utils.baseline import Baseline
            M = len(mets)
            r0 = fit_FSLModel(mrs, method='Newton', ppmlim=PPMLIM, baseline='poly, 2',
                              metab_groups=None, model='voigt')
            p0 = r0.params
            con, gam, sig, eps, ph0, ph1 = p0[:M], p0[M], p0[M + 1], p0[M + 2], p0[M + 3], p0[M + 4]
            nb = Baseline(mrs, PPMLIM, baseline, None).regressor.shape[1]
            x0 = np.r_[con, np.full(g, gam), np.full(g, sig), np.full(g, eps), ph0, ph1, np.zeros(nb)]
        res = fit_FSLModel(mrs, method='Newton', ppmlim=PPMLIM, baseline=baseline,
                           metab_groups=grp, model='voigt', x0=x0)
        c = np.asarray(res.getConc(scaling='raw', metab=mets), float)
        el = time.perf_counter() - t0
        if plot:
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                from fsl_mrs.utils.misc import FIDToSpec
                ppm = mrs.getAxes('ppmshift')
                y = FIDToSpec(mrs.FID)
                fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True)
                ax[0].plot(ppm, y.real, 'k', lw=0.8, label='data')
                ax[0].plot(ppm, np.real(FIDToSpec(res.pred)), 'r', lw=0.8, label='fit')
                ax[0].plot(ppm, np.real(FIDToSpec(res.baseline)), 'b', lw=0.8, label='baseline')
                ax[1].plot(ppm, y.real - np.real(FIDToSpec(res.pred)), 'k', lw=0.8)
                ax[0].legend(); ax[1].set_xlabel('ppm'); ax[0].set_xlim(10, -20)
                ax[0].set_title(f'{bname} basis, {baseline}'); ax[1].set_title('residual')
                fig.tight_layout(); fig.savefig(plot, dpi=110); plt.close(fig)
            except Exception as e:      # plotting must never stop the benchmark
                print(f'      (plot failed: {e!r})')
    if return_res:
        return c / mrs.fid_scaling, el, res, mrs
    return c / mrs.fid_scaling, el


def calibrate(cfg, outdir, try_conj=True):
    """Absolute scale K for cfg, or None if the noiseless amplitudes are not recovered."""
    bname, baseline, gmode = CONFIGS[cfg]
    F, mets, a = D['F_' + bname], D['mets_' + bname], D['aCal_' + bname]
    N = F.shape[0]
    t = np.arange(N) * D['dt']
    fid = (F @ a) * np.exp(-np.pi * D['LW_CAL'] * t)
    rng = np.random.default_rng(1)
    fid = fid + 1e-3 * np.abs(fid).max() * (rng.standard_normal(N) + 1j * rng.standard_normal(N))
    iP = mets.index('PCr')
    big = [i for i, n in enumerate(mets) if n in CAL_METS and a[i] > 0]
    variants = ((False, False), (True, True), (True, False), (False, True)) if try_conj else ((False, False),)
    png = os.path.join(outdir, 'calib_' + cfg.replace(' ', '_') + '.png')
    for iv, (cf, cb) in enumerate(variants):
        c, el = fit_one(fid, bname, baseline, gmode, cf, cb, plot=png if iv == 0 else None)
        rel_e = np.maximum(c[big] / c[iP], 1e-6)
        rel_t = a[big] / a[iP]
        err = float(np.sqrt(np.mean(np.log(rel_e / rel_t) ** 2)))
        print(f'  {cfg:20s} conj FID {int(cf)}, conj basis {int(cb)}: rms log-ratio error '
              f'{err:.4f} ({el:.1f} s)', flush=True)
        if iv == 0 and err >= 0.05:
            print('      metabolite   true/PCr   fitted/PCr')
            for i in range(len(mets)):
                if a[i] > 0:
                    print(f'      {mets[i]:10s} {a[i] / a[iP]:9.3f} {c[i] / c[iP]:11.3f}')
            print(f'      fit plot: {png}')
        if err < 0.05:
            if cf or cb:
                print(f'  {cfg}: only works with conjugation - export convention problem, skipped.')
                return None
            return a[iP] / c[iP]
    print(f'  {cfg}: calibration amplitudes not recovered - configuration SKIPPED.')
    return None


def work(task):
    ic, s, q, k, n, cfg, K = task
    bname, baseline, gmode = CONFIGS[cfg]
    try:
        c, el = fit_one(D['FID'][:, s, q, k, n], bname, baseline, gmode)
        amp = c * K                 # calibrated, in DSCALE units (PCr ~ 1)
        return task, amp, el, None
    except Exception as e:          # keep going, report at the end
        return task, None, np.nan, repr(e)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input', help='fsl_input.mat from fsl_export.m')
    ap.add_argument('--out', default=None, help='output .mat (default: fsl_results.mat next to input)')
    ap.add_argument('--configs', nargs='+', default=DEFAULT_CONFIGS, choices=list(CONFIGS))
    ap.add_argument('--scen', nargs='+', type=int, default=[4, 5, 6],
                    help='scenarios, 1-based as in MATLAB (default 4 5 6 = D E F)')
    ap.add_argument('--spline', default=None,
                    help="spline stiffness for the spline configs: very-stiff, stiff, moderate (default), "
                         "flexible, very-flexible, or an effective dimension, e.g. 10")
    ap.add_argument('--no-warm', action='store_true',
                    help='spline fits from the FSL-MRS default initialisation (no poly-2 warm start)')
    ap.add_argument('--no-cache', action='store_true',
                    help='recompute the spline penalty in every fit, as FSL-MRS does (for timing)')
    ap.add_argument('--resume', action='store_true',
                    help='skip the fits already saved in the output file (after an interrupted run)')
    ap.add_argument('--calib-only', action='store_true', help='run the calibration check only')
    ap.add_argument('--workers', type=int, default=max(1, (os.cpu_count() or 2) // 2))
    args = ap.parse_args()
    out = args.out or os.path.join(os.path.dirname(os.path.abspath(args.input)), 'fsl_results.mat')

    import fsl_mrs
    print(f'FSL-MRS {fsl_mrs.__version__}, {args.workers} workers')
    load_input(args.input)
    N, nS, nQ, nK, nN = D['FID'].shape
    REP = D['REP']
    cfgs = args.configs
    scen = [s - 1 for s in args.scen]

    global WARM_START, CACHE_LAMBDA
    WARM_START = not args.no_warm
    CACHE_LAMBDA = not args.no_cache
    if args.spline:
        for c in cfgs:
            if 'spline' in CONFIGS[c][1]:
                CONFIGS[c] = (CONFIGS[c][0], 'spline, ' + args.spline, CONFIGS[c][2])
    print('Calibration / convention check (near-noiseless):')
    outdir = os.path.dirname(os.path.abspath(out))
    K = {}
    for c in cfgs:     # conjugated variants only tried for the fast configurations
        K[c] = calibrate(c, outdir, try_conj='spline' not in CONFIGS[c][1])
    cfgs = [c for c in cfgs if K[c] is not None]
    if not cfgs:
        raise SystemExit('No configuration passed the calibration.')

    if args.calib_only:
        return
    DSCALE = float(np.squeeze(loadmat(args.input, variable_names=['DSCALE'])['DSCALE']))
    EST = np.full((nS, nQ, nK, len(cfgs), len(REP), nN), np.nan)
    T = np.full((nS, nQ, nK, len(cfgs), nN), np.nan)
    if args.resume and os.path.exists(out):      # keep the fits already saved in `out`
        old = loadmat(out)
        oldc = cellstr(old['cfgL'])
        E0, T0 = old['EST_L'], old['T_FIT']
        E0 = E0.reshape(E0.shape + (1,) * (6 - E0.ndim)); T0 = T0.reshape(T0.shape + (1,) * (5 - T0.ndim))
        for ic, c in enumerate(cfgs):
            if c in oldc:
                jo = oldc.index(c)
                EST[:, :, :, ic, :, :] = E0[:, :, :, jo, :, :]
                T[:, :, :, ic, :] = T0[:, :, :, jo, :]
        print(f'Resuming: {int(np.isfinite(T).sum())} fits already done in {out}')
    tasks = [(ic, s, q, k, n, c, K[c]) for ic, c in enumerate(cfgs)
             for s in scen for q in range(nQ) for k in range(nK) for n in range(nN)
             if not np.isfinite(T[s, q, k, ic, n])]
    # slow (spline) tasks first, so the pool finishes evenly
    tasks.sort(key=lambda x: 'spline' not in CONFIGS[x[5]][1])
    errors = []

    def save():
        savemat(out, {'EST_L': EST, 'cfgL': np.array(cfgs, dtype=object), 'T_FIT': T,
                      'K_CAL': np.array([K[c] for c in cfgs]), 'SCEN_RUN': np.array(args.scen),
                      'CFG_BASELINE': np.array([CONFIGS[c][1] for c in cfgs], dtype=object)})

    t_start = time.time()
    with Pool(args.workers) as pool:
        for i, (task, amp, el, err) in enumerate(pool.imap_unordered(work, tasks), 1):
            ic, s, q, k, n, cfg, _ = task
            if err is None:
                W = D['W_' + CONFIGS[cfg][0]]
                EST[s, q, k, ic, :, n] = W @ (amp / DSCALE)
                T[s, q, k, ic, n] = el
            else:
                errors.append((cfg, s + 1, q + 1, k + 1, err))
            if i % 50 == 0 or i == len(tasks):
                print(f'  {i}/{len(tasks)} fits, {time.time() - t_start:.0f} s', flush=True)
                save()
    save()
    print(f'Saved {out}')
    for ic, c in enumerate(cfgs):
        tt = T[:, :, :, ic, :]
        print(f'  {c:20s} median {np.nanmedian(tt):6.2f} s/fit (IQR {np.nanpercentile(tt, 25):.2f}-'
              f'{np.nanpercentile(tt, 75):.2f}) with {args.workers} parallel workers')
    if errors:
        print(f'{len(errors)} fits failed, first: {errors[0]}')


if __name__ == '__main__':
    main()
