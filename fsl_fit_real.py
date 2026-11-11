#!/usr/bin/env python
"""FSL-MRS fit of a real 31P spectrum exported by fit_real_3way.m.

    $FSLDIR/bin/python fsl_fit_real.py /mnt/c/.../real_fits/fsl_real_input.mat

Uses the configurations and fitting code of fsl_benchmark.py (same folder):
calibration check on a near-noiseless copy of the in-house fit, then the fit,
then the curves for the figures (data, fit, baseline and one curve per
metabolite, fitted phase removed). Writes fsl_real_results.mat next to the input.
"""
import os
for _v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
    os.environ.setdefault(_v, '1')          # single thread, set before numpy is imported
import argparse  # noqa: E402
import sys  # noqa: E402

import numpy as np  # noqa: E402
from scipy.io import savemat  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fsl_benchmark as fb          # noqa: E402

DEFAULT = ['FSL-NEW spline grpA']   # add 'FSL-NEW spline' 'FSL-OLD spline' with --configs if wanted


def curves(res, mrs, mets, g):
    from fsl_mrs.utils.misc import FIDToSpec
    M = len(mets)
    p = res.params
    phi0, phi1 = p[M + 3 * g], p[M + 3 * g + 1]
    ph = np.exp(1j * (phi0 + phi1 * np.ravel(mrs.frequencyAxis)))   # undo the fitted phase for display
    fs = mrs.fid_scaling
    sp = lambda fid: FIDToSpec(np.ravel(np.asarray(fid)).copy()) * ph / fs   # noqa: E731
    comps = np.column_stack([sp(res.predictedFID(mrs, mode=m, noBaseline=True)) for m in mets])
    return sp(mrs.FID), sp(res.pred), sp(res.baseline), comps


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('input')
    ap.add_argument('--configs', nargs='+', default=DEFAULT, choices=list(fb.CONFIGS))
    ap.add_argument('--repeat', type=int, default=1,
                    help='fits repeated for the timing comparison (time_s = median; single thread)')
    args = ap.parse_args()
    out = os.path.join(os.path.dirname(os.path.abspath(args.input)), 'fsl_real_results.mat')

    fb.load_input(args.input)
    fid = fb.D['FID'][:, 0, 0, 0, 0]
    outdir = os.path.dirname(out)
    R = []
    for cfg in args.configs:
        bname, baseline, gmode = fb.CONFIGS[cfg]
        print(f'{cfg}: calibration ...', flush=True)
        K = fb.calibrate(cfg, outdir, try_conj=False)
        if K is None:
            print(f'{cfg}: skipped (calibration failed)')
            continue
        times = []
        for _ in range(max(1, args.repeat)):     # first fit also computes the cached spline penalty
            c, el, res, mrs = fb.fit_one(fid, bname, baseline, gmode, return_res=True)
            times.append(el)
        el = float(np.median(times))
        mets = fb.D['mets_' + bname]
        grp = fb.groups_for(mets, gmode)
        g = 1 if grp is None else max(grp) + 1
        amp = c * K
        data, fit, base, comps = curves(res, mrs, mets, g)
        try:
            crlb = np.asarray(res.getUncertainties(type='percentage', metab=mets), float)
        except Exception:                     # older/newer API differences
            crlb = np.full(len(mets), np.nan)
        R.append({'name': cfg, 'basis': bname, 'mets': np.array(mets, dtype=object), 'amp': amp,
                  'est': fb.D['W_' + bname] @ amp, 'crlb_pct': crlb, 'time_s': el,
                  'time_all': np.array(times),
                  'data': data, 'fit': fit, 'base': base, 'comps': comps})
        print(f'{cfg}: done, median {el:.1f} s over {len(times)} fit(s): {np.round(times, 2)}', flush=True)
    cell = np.empty((len(R),), dtype=object)
    for i, r in enumerate(R):
        cell[i] = r
    savemat(out, {'R': cell, 'REP': np.array(fb.D['REP'], dtype=object)})
    print(f'Saved {out}')


if __name__ == '__main__':
    main()
