%% FIT_REAL_3WAY  (4-way: in-house, LCModel, AMARES, FSL-MRS)
% Fit one real 31P spectrum with the in-house fitter, LCModel, AMARES (OXSA) and
% FSL-MRS, with the corrected (Ren) basis (optionally also the original basis),
% and draw the fits.
%
% Run sections 0-4 in MATLAB. Section 4 writes the FSL-MRS input and prints
% the WSL command; run it, then run section 5 (figures + table), which also
% works without FSL-MRS results (it then shows two fitters).
%
% Needs in the workspace: basis_real (original basis struct) and `spectrum`
% (complex, 1024 points, from your TWIX code), or specL + VOXELS for an aligned
% multi-voxel average. Files: fit_31P_spectrum_v4_cs,
% make_basis_ren2020, and fsl_benchmark.py + fsl_fit_real.py in the same folder.

%% 0. settings -------------------------------------------------------------
DWELL_ORIG = 1e-4;     % s, dwell time of the ACQUIRED FID. You keep 1024 of 2048 spectral
                       % points, so the kept spectrum has dwell 2*DWELL_ORIG; it must equal
                       % the basis dwell (2e-4 s). Check it in the TWIX header.
FIT_OLD    = false;    % original basis on real data (off: the basis comparison is done on simulations)
% OXSA_DIR   = '';       % OXSA toolbox folder (contains main/ and utils/) for AMARES; '' = skip AMARES
AMARES_PHASE = 'nad_aatp'; % own phase per peak, NAD/NADH/UDP tied to alpha-ATP ('all' = one phase)
AMARES_SKIP  = 2;      % first FID points left out by AMARES (0 = none, as usual for AMARES)
AMARES_BROAD = [];     % ppm: optional broad peaks for a baseline hump ([] = none, as usual for AMARES)
AMARES_LINESHAPE = 'voigt';  % as in the benchmark: Lorentz x Gauss with one shared Gaussian width
VOXELS     = [165 157 164];       % e.g. []: aligned average of these columns of specL
                       % (each voxel's PCr frequency and phase aligned before summing);
                       % [] = fit `spectrum` as given
OUTDIR     = fullfile(pwd, 'real_fits_S2_avg4_pref/');
LCM_BIN    = 'LCModel_exe\LCModel.exe';
LCM_DIR    = ''; KEY = []; OWNER = '';
ECHOT      = 0.35;     % ms
LCM_OPTS = { ...
    'PPMST', 10, 'PPMEND', -20, 'NSIMUL', 0, 'NRATIO', 0, 'DKNTMN', 2, 'RFWHM', 3, ...
    'SDDEGZ', 2, 'SDDEGP', 0.5, 'DEGZER', 0, 'DEGPPM', 0, 'DESDT2', 10, 'DEEXT2', 10, ...
    'DESDSH', 0.01, 'DOECC', 'F', 'DOWS', 'F', 'DOREFS', 'F', 'NCALIB', 0, 'NEACH', 999, ...
    'LCOORD', 9, 'LTABLE', 7, 'NAMREL', 'PCr', 'CONREL', 1};
CHUSE1 = {'PCr','PiIn','aATP','gATP'};
N_TIMING      = 1;     % fits repeated per fitter for the timing comparison (e.g. 5); time_s = median
SINGLE_THREAD = false; % true: maxNumCompThreads(1), so the fitters are timed on one core
                       % (FSL-MRS runs single-threaded already; LCModel is single-threaded)
CHCOMB = {'NADH+NADp'};

if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end
if SINGLE_THREAD, nthr0 = maxNumCompThreads(1); end   % restored at the end of section 3b
bNew = make_basis_ren2020(basis_real, 'normalise_all', true);
bOld = basis_real;
B0 = bNew.B0_MHz; t = bNew.t(:); t = t - t(1); N = numel(t); dt = bNew.t(2) - bNew.t(1);
ppm = bNew.ppm(:);
REP = {'PCr','Pi_intra','PE','PC','GPC','GPE','gATP','aATP','NADplus','NADH','tNAD','UDP'};

%% 1. preprocessing: orientation, PCr at 0 ppm, phase -----------------------
if isempty(VOXELS)
    spec = double(spectrum(:));
else
    spec = align_average(double(specL(:, VOXELS)), ppm, B0, t);
end
if numel(spec) ~= N, error('spectrum has %d points, the basis %d.', numel(spec), N); end
if abs(2*DWELL_ORIG - dt) > 1e-9
    error(['Kept spectrum dwell %.3g s does not match the basis dwell %.3g s. ' ...
           'Check DWELL_ORIG or crop/zero-fill differently.'], 2*DWELL_ORIG, dt);
end
dppm = ppm(2) - ppm(1);
fid = ifft(ifftshift(spec));

% orientation: beta-ATP (-16.1 ppm) and alpha-ATP (-7.55) must be upfield of PCr
p = pcr_index(fftshift(fft(fid)));
sc_dn = win_mag(fftshift(fft(fid)), p, [-16.1 -7.55], dppm);
sc_up = win_mag(fftshift(fft(fid)), p, [16.1 7.55], dppm);
if sc_up > sc_dn
    fid = conj(fid);
    fprintf('Spectrum was frequency-reversed: FID conjugated.\n');
end

% PCr to 0 ppm (sub-bin)
S0 = fftshift(fft(fid));
p = pcr_index(S0);
ppm_p = ppm(1) + (p - 1)*dppm;
fid = fid .* exp(-1i*2*pi*ppm_p*B0*t);
S1 = fftshift(fft(fid));
fprintf('PCr found at %.3f ppm on the input axis; moved to 0.\n', ppm_p);

% phase: zero order + linear (t0), from the coherence of the main peak integrals
[t0e, ph0e] = estimate_phase(S1, ppm, B0);
S2 = S1 .* exp(-1i*(ph0e + 2*pi*B0*t0e*ppm));
fprintf('Phase estimate: phi0 = %.1f deg, t0 = %.0f us\n', ph0e*180/pi, t0e*1e6);

% check peak positions (tests the bandwidth assumption)
CHK = {'gATP', -2.48, [-3.2 -1.8]; 'aATP', -7.55, [-8.2 -6.9]; 'Pi', 4.9, [4.2 5.6]; 'bATP', -16.15, [-17 -15.3]};
for k = 1:size(CHK,1)
    m = ppm >= CHK{k,3}(1) & ppm <= CHK{k,3}(2);
    [~, j] = max(real(S2) .* m);
    fprintf('  %-5s peak at %6.2f ppm (expected ~%6.2f)\n', CHK{k,1}, ppm(j), CHK{k,2});
    if abs(ppm(j) - CHK{k,2}) > 0.3
        warning(['%s is %.2f ppm off. If gamma/alpha-ATP and Pi are right, %s is probably just ' ...
                 'weak (e.g. outside the excitation bandwidth).'], CHK{k,1}, ppm(j) - CHK{k,2}, CHK{k,1});
    end
end

DSCALE = 100 / max(real(S2(abs(ppm) < 0.5)));
spec_p = S2 * DSCALE;                       % phased, referenced, scaled spectrum
fid_p  = ifft(ifftshift(spec_p));
figure('Name','Preprocessed spectrum','Color','w');
plot(ppm, real(spec_p), 'k'); set(gca,'XDir','reverse'); xlim([-20 10]); grid on;
xlabel('ppm'); title('Real part after referencing and phasing');

%% 2. in-house fitter ----------------------------------------------------------
LGW = {'gATP', {'gATP'}, 15; 'aATP', {'aATP'}, 15; 'Pi', {'Pi_intra','Pi_extra'}, 15};
BAS = struct('name', {'NEW','OLD'}, 'b', {bNew, bOld}, 'nadprior', {[0 5], [15 10]});
if ~FIT_OLD, BAS = BAS(1); end
FITS = struct('fitter', {}, 'basis', {}, 'ppm', {}, 'data', {}, 'fit', {}, 'bg', {}, ...
              'comp', {}, 'compnames', {}, 'est', {}, 'crlb', {}, 'time_s', {}, 'time_all', {});
for i = 1:numel(BAS)
    b = BAS(i).b; M = numel(b.metabolites);
    AB = [-0.1*ones(1,M); Inf(1,M)]; AB(:, strcmp(b.metabolites,'MP_broad')) = [0; 1e-12];
    nn = b.metabolites(ismember(b.metabolites, {'NADH','NADplus','UDPG'}));
    P = NaN(4,2); P(4,:) = BAS(i).nadprior;
    args = {'ppm_range',[-20 10], 'plot',false, 't0',0, 't0_prior_sd',25e-6, ...
            'atp_shift_prior_hz',[0 8; 0 8], 'nad_rel_aatp',true, 'lw_groups',LGW, ...
            'amp_bounds',AB, 'nonneg_mets',nn, 'fsh_grp_prior_hz',P, ...
            'n_baseline_knots',20, 'baseline_curv_sd',1, 'fit_eta',true};
    tt = zeros(N_TIMING, 1);
    for kk = 1:N_TIMING, tic; r = fit_31P_spectrum_v4_cs(spec_p, b, args{:}); tt(kk) = toc; end
    el = median(tt);
    F = struct(); F.fitter = 'In-house'; F.basis = BAS(i).name; F.ppm = r.ppm(:);
    F.data = r.spec_obs(:); F.fit = r.spec_fit(:); F.bg = r.baseline_full(:);
    F.comp = r.spec_each; F.compnames = r.metabolites(:)';
    F.est = report_amounts(r.amplitudes_raw(:), b, r.metabolites);
    F.crlb = containers.Map(r.metabolites, num2cell(r.metrics.CRLB_pct(:)'));
    F.time_s = el; F.time_all = tt;
    FITS(end+1) = F; %#ok<SAGROW>
    fprintf('In-house, %s basis: t0 %.0f us, eta %.2f, NAD offset vs aATP %+.1f Hz, %.1f s\n', ...
        BAS(i).name, r.t0*1e6, r.eta, r.nad_rel_aatp_hz, el);
    if i == 1, rNew = r; end
end
save(fullfile(OUTDIR, 'real_fits.mat'), 'FITS', 'spec_p', 'ppm', 'DSCALE', 't0e', 'ph0e', 'REP');

%% 3. LCModel -------------------------------------------------------------------
if isempty(LCM_BIN)
    warning('LCM_BIN not set: LCModel skipped.');
else
    MAP = {'PCr','PCr'; 'Pi_intra','PiIn'; 'Pi_extra','PiEx'; 'PE','PE'; 'PC','PC'; ...
           'GPC','GPC'; 'GPE','GPE'; 'gATP','gATP'; 'aATP','aATP'; 'bATP','bATP'; ...
           'NADH','NADH'; 'NADplus','NADp'; 'UDPG','UDPG'; 'MP_broad','MPbr'};
    short = @(m) MAP{strcmp(MAP(:,1), m), 2};
    pfit = polyfit(((0:N-1)' - floor(N/2))/(N*dt)/B0, ppm, 1); PPMCEN = pfit(2);
    fwhmba = mean(bNew.prior.lw_T2_hz) / B0;
    fraw = fullfile(OUTDIR, 'real.RAW');
    write_raw(fraw, fid_p, true, ECHOT, B0);        % convention verified in lcm_benchmark (RAW conj, basis conv 3)
    for i = 1:numel(BAS)
        b = BAS(i).b;
        [Fb, mets] = basis_fids(b);
        keep = ~strcmp(mets, 'MP_broad'); Fb = Fb(:, keep); mets = mets(keep);
        sn = cellfun(short, mets, 'UniformOutput', false);
        fbas = fullfile(OUTDIR, sprintf('%s.BASIS', BAS(i).name));
        write_basis(fbas, Fb, sn, dt, B0, ECHOT, fwhmba, 3, PPMCEN);
        stem = fullfile(OUTDIR, sprintf('lcm_%s', BAS(i).name));
        tt = zeros(N_TIMING, 1);       % includes writing the control file and starting LCModel, as per voxel in CSI
        for kk = 1:N_TIMING
            tic; [tab, sd] = run_lcm(stem, fraw, fbas, dt, B0, N, PPMCEN, ECHOT, LCM_OPTS, CHUSE1, CHCOMB, ...
                                     KEY, OWNER, LCM_BIN, LCM_DIR); tt(kk) = toc;
        end
        el = median(tt);
        C = parse_coord([stem '.coord']);
        amp = zeros(numel(b.metabolites), 1); crl = containers.Map();
        for m = 1:numel(mets)
            if isKey(tab, sn{m}), amp(strcmp(b.metabolites, mets{m})) = tab(sn{m}); crl(mets{m}) = sd(sn{m}); end
        end
        F = struct(); F.fitter = 'LCModel'; F.basis = BAS(i).name; F.ppm = C.ppm(:);
        F.data = C.data(:); F.fit = C.fit(:); F.bg = C.bg(:);
        F.comp = C.curves - C.bg(:); F.compnames = cellfun(@(s) long_name(s, MAP), C.mets, 'UniformOutput', false);
        F.est = report_amounts(amp, b, b.metabolites); F.crlb = crl; F.time_s = el; F.time_all = tt;
        FITS(end+1) = F; %#ok<SAGROW>
        fprintf('LCModel, %s basis: %d metabolites in the table, %.1f s\n', BAS(i).name, tab.Count, el);
    end
    save(fullfile(OUTDIR, 'real_fits.mat'), 'FITS', '-append');
end

%% 3b. AMARES (OXSA) -------------------------------------------------------------
% Prior knowledge generated from the same corrected basis, with the four
% UDP-sugars as separate peaks (as in Ren et al.) and NAD+ as two doublets
% (exact AB positions); see make_amares_pk.
% if isempty(OXSA_DIR)
%     warning('OXSA_DIR not set: AMARES skipped.');
% else
%     addpath(genpath(fullfile(OXSA_DIR, 'main'))); addpath(genpath(fullfile(OXSA_DIR, 'utils')));
    bAm = make_basis_ren2020(basis_real, 'normalise_all', true, 'udp_mode', 'separate');
    [pkA, infoA] = make_amares_pk(bAm, 'nad_quartet', 'two_doublets', 'phase_mode', AMARES_PHASE, 'broad_ppm', AMARES_BROAD, 'lineshape', AMARES_LINESHAPE, 'write', fullfile(OUTDIR, 'PK_7T_Brain_Ren2020_31P.m'));
    % AMARES gets the data WITHOUT the first-order phase correction, and models the
    % acquisition delay itself (beginTime = t0e): its lines are evaluated from
    % t0e on, like the acquired FID. A first-order phase applied to the spectrum
    % is a circular shift of the FID, which creates the broad baseline roll that
    % a model without baseline cannot fit. Only zero-order phase and scaling here.
    spec_a = S1 .* exp(-1i*ph0e) * DSCALE;
    tt = zeros(N_TIMING, 1);
    for kk = 1:N_TIMING
        tic; A = amares_fit_spectrum(spec_a, bAm, pkA, infoA, 'skip_points', AMARES_SKIP, 'beginTime', t0e); tt(kk) = toc;
    end
    el = median(tt);
    lin = exp(-1i*2*pi*B0*t0e*ppm);             % curves to the frame of spec_p (same rotation for data and fit)
    A.data = A.data .* lin; A.fit = A.fit .* lin; A.bg = A.bg .* lin; A.comp = A.comp .* lin;
    [est, comp, cn, crl] = amares_to_rep(A, REP);
    F = struct(); F.fitter = 'AMARES'; F.basis = 'NEW'; F.ppm = ppm;
    F.data = real(A.data); F.fit = real(A.fit);   % fit = all AMARES peaks (+ broad peaks / skipped points if used)
    F.bg = real(A.bg);                             % zero unless AMARES_BROAD / AMARES_SKIP are used
    F.comp = real(comp); F.compnames = cn; F.est = est; F.crlb = crl; F.time_s = el; F.time_all = tt;
    FITS(end+1) = F;
    save(fullfile(OUTDIR, 'real_fits.mat'), 'FITS', '-append');
    fprintf('AMARES: %d peaks, %.1f s\n', numel(infoA.rows), el);
% end
if SINGLE_THREAD, maxNumCompThreads(nthr0); end

%% 4. FSL-MRS input (fit in WSL) ----------------------------------------------
out = struct('FID', reshape(fid_p, [N 1 1 1 1]), 'dt', dt, 'B0', B0, 'DSCALE', 1, 'LW_CAL', 8, ...
             'dims', [N 1 1 1 1]);
out.REP = REP;
aN = max(rNew.amplitudes_raw(:), 0); aN(strcmp(bNew.metabolites, 'MP_broad')) = 0;
bs = {'OLD', bOld; 'NEW', bNew};
for i = 1:2
    b = bs{i,2};
    [Fb, mets] = basis_fids(b);
    keep = ~strcmp(mets, 'MP_broad'); Fb = Fb(:, keep); mets = mets(keep);
    aCal = zeros(numel(mets), 1);
    for m = 1:numel(mets), j = strcmp(bNew.metabolites, mets{m}); if any(j), aCal(m) = aN(j); end, end
    out.(['F_' bs{i,1}]) = Fb; out.(['mets_' bs{i,1}]) = mets;
    out.(['W_' bs{i,1}]) = rep_matrix(b, mets, REP); out.(['aCal_' bs{i,1}]) = aCal;
end
fin = fullfile(OUTDIR, 'fsl_real_input.mat');
save(fin, '-struct', 'out', '-v7');
fprintf('\nNow in WSL (folder with fsl_benchmark.py and fsl_fit_real.py):\n  $FSLDIR/bin/python fsl_fit_real.py %s --repeat %d\nthen run section 5.\n', wsl_path(fin), N_TIMING);

%% 5. figures and table --------------------------------------------------------
L = load(fullfile(OUTDIR, 'real_fits.mat'));
FT = L.FITS; ppmA = L.ppm;
ffsl = fullfile(OUTDIR, 'fsl_real_results.mat');
if exist(ffsl, 'file')
    G = load(ffsl); Rc = G.R; if ~iscell(Rc), Rc = num2cell(Rc); end
    for i = 1:numel(Rc)
        r = Rc{i};
        F = FT(1); F.fitter = ['FSL-MRS' strrep(strrep(char(r.name), 'FSL-NEW', ''), 'FSL-OLD', '')];
        F.basis = char(r.basis); F.ppm = ppmA;
        F.data = real(r.data(:)); F.fit = real(r.fit(:)); F.bg = real(r.base(:));
        F.comp = real(r.comps); F.compnames = cellstr(r.mets); F.compnames = F.compnames(:)';
        F.est = r.est(:)'; F.crlb = containers.Map(F.compnames, num2cell(r.crlb_pct(:)'));
        F.time_s = r.time_s;
        if isfield(r, 'time_all'), F.time_all = r.time_all(:); else, F.time_all = r.time_s; end
        FT(end+1) = F; %#ok<SAGROW>
    end
else
    fprintf('(no FSL-MRS results yet: %s)\n', ffsl);
end

% table: per molecule relative to PCr
iP = strcmp(REP, 'PCr'); iNp = strcmp(REP, 'NADplus'); iNh = strcmp(REP, 'NADH');
fprintf('\n%-22s %-5s | %8s %8s %8s %8s %8s %8s | %7s %7s | %6s\n', 'fitter', 'basis', 'NAD+/PCr', 'NADH/PCr', ...
    'tNAD/PCr', 'UDP/PCr', 'aATP/PCr', 'Pi/PCr', 'CRLB N+', 'CRLB NH', 'ratio');
for i = 1:numel(FT)
    e = FT(i).est / FT(i).est(iP);
    fprintf('%-22s %-5s | %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f | %6.0f%% %6.0f%% | %6.1f\n', FT(i).fitter, FT(i).basis, ...
        e(iNp), e(iNh), e(strcmp(REP,'tNAD')), e(strcmp(REP,'UDP')), e(strcmp(REP,'aATP')), e(strcmp(REP,'Pi_intra')), ...
        getc(FT(i).crlb, 'NADplus'), getc(FT(i).crlb, 'NADH'), safe_ratio(e(iNp), e(iNh)));
end

% All fitters on the same ppm grid and scale (curves interpolated to the grid
% and scaled by regressing each fitter's copy of the data onto spec_p). Each
% fit is drawn over the data AS THAT FITTER PHASED IT: LCModel and FSL-MRS
% remove their own fitted phase from the data, so drawing their fits over
% spec_p would show dispersion-shaped residuals that are not in their fits.
% Baselines are drawn, not subtracted.
yref = real(L.spec_p(:));
for i = 1:numel(FT), FT(i) = to_common(FT(i), ppmA, yref); end
nrm = max(yref(abs(ppmA) < 0.5));
ZOOM = [-10.5 -5.5];
SHOW = {'aATP','NADplus','NADH','UDPG'}; COL = [0.55 0.55 0.55; 0.85 0.33 0.10; 0 0.45 0.74; 0.47 0.67 0.19];

% Figure 1: corrected basis, one column per fitter: full spectrum + NAD region
sel = find(strcmp({FT.basis}, 'NEW'));
figure('Name','Real spectrum: fits, corrected basis','Color','w','Position',[30 30 420*numel(sel) 720]);
for k = 1:numel(sel)
    F = FT(sel(k)); x = ppmA;
    subplot(2, numel(sel), k); hold on;
    plot(x, F.data/nrm, 'Color', [0.6 0.6 0.6]); plot(x, F.fit/nrm, 'r', 'LineWidth', 1.1);
    plot(x, F.bg/nrm, 'b', 'LineWidth', 1.1);
    plot(x, (F.data - F.fit)/nrm - 0.25, 'k');
    set(gca,'XDir','reverse'); xlim([-20 10]); grid on; box off; xlabel('ppm');
    title(sprintf('%s (%s basis)', F.fitter, F.basis));
    if k == 1, legend({'data','fit','baseline','residual (offset)'}, 'Location','northwest','FontSize',7); end
    subplot(2, numel(sel), numel(sel) + k);
    nad_panel(F, x, nrm, ZOOM, SHOW, COL, iP, iNp, iNh, '');
    if k == 1, legend([{'data','fit','baseline'}, SHOW], 'Location','northwest','FontSize',7); end
end
sgtitle('Real spectrum, corrected basis (common grid and scale; data as phased by each fitter)', 'FontWeight','bold');

% Figure 2: original vs corrected basis in the NAD region, per fitter
fitters = unique({FT.fitter}, 'stable');
if any(strcmp({FT.basis}, 'OLD'))
    figure('Name','Real spectrum: original vs corrected basis','Color','w','Position',[60 60 380*numel(fitters) 640]);
    BN = {'OLD','NEW'};
    for k = 1:numel(fitters)
        for bb = 1:2
            j = find(strcmp({FT.fitter}, fitters{k}) & strcmp({FT.basis}, BN{bb}), 1);
            subplot(2, numel(fitters), (bb-1)*numel(fitters) + k);
            if isempty(j), axis off; continue; end
            nad_panel(FT(j), ppmA, nrm, ZOOM, SHOW, COL, iP, iNp, iNh, sprintf('%s, %s: ', FT(j).fitter, lower(BN{bb})));
        end
    end
    sgtitle('NAD region: original (top) vs corrected (bottom) basis', 'FontWeight','bold');
end

% Figure 3: summary bars
figure('Name','Real spectrum: NAD estimates','Color','w','Position',[90 90 1100 380]);
lab = strcat({FT.fitter}, {' / '}, lower({FT.basis}));
Q = {'NADplus','NADH','tNAD','UDP'};
V = nan(numel(FT), numel(Q));
for i = 1:numel(FT), for q = 1:numel(Q), V(i,q) = FT(i).est(strcmp(REP,Q{q})) / FT(i).est(iP); end, end
subplot(1,2,1); bar(V); set(gca,'XTickLabel', lab, 'XTickLabelRotation', 30); grid on; box off;
ylabel('per molecule / PCr'); legend({'NAD+','NADH','tNAD','UDP(G)'}, 'Location','best');
subplot(1,2,2); bar(arrayfun(@(i) safe_ratio(V(i,1), V(i,2)), 1:numel(FT)));
set(gca,'XTickLabel', lab, 'XTickLabelRotation', 30, 'YScale','log'); grid on; box off;
ylabel('NAD+/NADH (log)'); yline(5, 'k:', 'healthy > ~5');

%% ===========================================================================
function [est, comp, cn, crl] = amares_to_rep(A, REP)
% AMARES per-molecule amounts on the REP list; UDP-sugars summed to UDP(G).
isu = ~cellfun(@isempty, regexp(A.mets, '^UDP'));
g = @(n) sum(A.per_mol(strcmp(A.mets, n)));
est = nan(1, numel(REP));
for j = 1:numel(REP)
    switch REP{j}
        case 'tNAD', est(j) = g('NADH') + g('NADplus');
        case 'UDP',  est(j) = sum(A.per_mol(isu));
        otherwise
            if any(strcmp(A.mets, REP{j})), est(j) = g(REP{j}); end
    end
end
keep = ~isu;
comp = [A.comp(:, keep), sum(A.comp(:, isu), 2)];
cn = [A.mets(keep), {'UDPG'}];
crl = containers.Map(A.mets, num2cell(A.crlb_pct));
end

function spec = align_average(SV, ppm, B0, t)
% Aligned average of voxel spectra (columns of SV, fftshifted spectra): each
% voxel's PCr is moved to the frequency of the first voxel's PCr and its PCr
% peak rotated to the first voxel's phase, then all are averaged.
dppm = ppm(2) - ppm(1);
nv = size(SV, 2);
F = zeros(size(SV));
for v = 1:nv
    fid = ifft(ifftshift(SV(:, v)));
    p = pcr_index(fftshift(fft(fid)));
    if v == 1, p1 = p; end
    fid = fid .* exp(-1i*2*pi*(p - p1)*dppm*B0*t);          % frequency to voxel 1
    S = fftshift(fft(fid));
    ph = pcr_phase(S, round(p1));
    if v == 1, ph1 = ph; end
    F(:, v) = fid * exp(-1i*(ph - ph1));                     % zero-order phase to voxel 1
    fprintf('  voxel %d: PCr shift %+6.1f Hz, phase %+6.1f deg\n', v, (p - p1)*dppm*B0, (ph - ph1)*180/pi);
end
spec = fftshift(fft(mean(F, 2)));
naive = fftshift(fft(mean(cell2mat(arrayfun(@(v) ifft(ifftshift(SV(:,v))), 1:nv, 'UniformOutput', false)), 2)));
q = @(S) max(abs(S)) / noise_sd(S, ppm);
fprintf('  PCr peak/noise: single voxel %.0f, plain mean %.0f, aligned mean %.0f (%d voxels)\n', ...
    q(SV(:,1)), q(naive), q(spec), nv);
end

function ph = pcr_phase(S, ip)
w = max(1, ip-3):min(numel(S), ip+3);
ph = angle(sum(S(w)));
end

function sd = noise_sd(S, ppm)
r = [diff(S(ppm > 8.5)); diff(S(ppm < -18.5))];
sd = std([real(r); imag(r)]) / sqrt(2);
end

function F = to_common(F, x, yref)
% Fitter curves onto the common grid x, scaled to the reference data.
xs = F.ppm(:);
if numel(xs) ~= numel(x) || any(abs(xs - x) > 1e-9)
    [xs, o] = sort(xs);
    g = @(y) interp1(xs, y(o, :), x, 'linear', NaN);
    F.data = g(F.data(:)); F.fit = g(F.fit(:)); F.bg = g(F.bg(:));
    if ~isempty(F.comp), F.comp = g(F.comp); end
    F.ppm = x;
end
m = ~isnan(F.data) & abs(x) < 12;
k = (F.data(m)' * yref(m)) / (F.data(m)' * F.data(m));
F.data = k*F.data; F.fit = k*F.fit; F.bg = k*F.bg; F.comp = k*F.comp;
end

function nad_panel(F, x, nrm, ZOOM, SHOW, COL, iP, iNp, iNh, pre)
hold on;
m = x >= ZOOM(1) & x <= ZOOM(2);
plot(x(m), F.data(m)/nrm, 'Color', [0.6 0.6 0.6]);
plot(x(m), F.fit(m)/nrm, 'r', 'LineWidth', 1.1);
plot(x(m), F.bg(m)/nrm, 'b', 'LineWidth', 1.1);
for c = 1:numel(SHOW)
    j = find(strcmp(F.compnames, SHOW{c}), 1);
    if ~isempty(j), plot(x(m), (F.bg(m) + F.comp(m, j))/nrm, 'Color', COL(c,:), 'LineWidth', 1.3); end
end
e = F.est / F.est(iP);
set(gca,'XDir','reverse'); xlim(ZOOM); grid on; box off; xlabel('ppm');
title(sprintf('%sNAD+/PCr %.3f, NADH/PCr %.3f, ratio %.1f', pre, e(iNp), e(iNh), safe_ratio(e(iNp), e(iNh))), 'FontSize', 8);
end

function p = pcr_index(S)
a = abs(S); a = conv(a, ones(5,1)/5, 'same');
[~, p] = max(a);
if p > 1 && p < numel(a)               % parabolic sub-bin refinement
    y = a(p-1:p+1); p = p + 0.5*(y(1) - y(3)) / (y(1) - 2*y(2) + y(3));
end
end

function v = win_mag(S, p, cen, dppm)
v = 0;
for c = cen
    i0 = round(p + c/dppm); w = round(0.6/dppm);
    i = max(1, i0 - w):min(numel(S), i0 + w);
    if ~isempty(i), v = v + sum(abs(S(i))); end
end
end

function [t0, ph0] = estimate_phase(S, ppm, B0)
% PCr, gamma-ATP, alpha-ATP, Pi, PE. Not beta-ATP: it may be missing when the
% excitation bandwidth does not reach it, and a noise window would spoil the estimate.
W = {[-0.6 0.6], [-3.2 -1.8], [-8.2 -6.9], [4.3 5.6], [6.4 7.2]};
msk = cellfun(@(w) ppm >= w(1) & ppm <= w(2), W, 'UniformOutput', false);
J = @(t0) -abs(sum(cellfun(@(m) unit(sum(S(m) .* exp(-1i*2*pi*B0*t0*ppm(m)))), msk)));
tg = (-1500:2:1500)*1e-6; Jg = arrayfun(J, tg);
[~, k] = min(Jg);
t0 = fminbnd(J, tg(max(k-2,1)), tg(min(k+2,end)));
z = sum(cellfun(@(m) unit(sum(S(m) .* exp(-1i*2*pi*B0*t0*ppm(m)))), msk));
ph0 = angle(z);
end

function u = unit(z), u = z / max(abs(z), eps); end

function [nrm, x] = norm_pcr(F)
x = F.ppm(:);
nrm = max(F.fit(abs(x) < 0.5));
if isempty(nrm) || nrm <= 0, nrm = max(abs(F.fit)); end
end

function c = getc(M, k)
if isKey(M, k), c = M(k); else, c = NaN; end
end

function r = safe_ratio(a, b)
if b <= 0, r = Inf; else, r = a / b; end
end

function n = long_name(s, MAP)
j = strcmp(MAP(:,2), s);
if any(j), n = MAP{find(j,1), 1}; else, n = s; end
end

function [F, mets] = basis_fids(b)
t = b.t(:); B0 = b.B0_MHz; mets = b.metabolites(:)';
F = zeros(numel(t), numel(mets));
for m = 1:numel(mets)
    pk = b.prior(strcmp(b.prior.metabolite, mets{m}), :);
    for r = 1:height(pk)
        F(:, m) = F(:, m) + pk.rel_amp(r) * exp(1i*2*pi*pk.cs_ppm(r)*B0*t - pi*pk.lw_T2_hz(r)*t);
    end
end
end

function write_basis(fname, F, names, dt, B0, echot, fwhmba, conv, ppmcen)
fid = fopen(fname, 'w');
fprintf(fid, ' $SEQPAR\n FWHMBA = %.5f,\n HZPPPM = %.5f,\n ECHOT = %.2f,\n SEQ = ''PULSEACQ''\n $END\n', ...
        fwhmba, B0, echot);
fprintf(fid, ' $BASIS1\n IDBASI = ''real'',\n FMTBAS = ''(2E15.6)'',\n BADELT = %.6e,\n NDATAB = %d\n $END\n', ...
        dt, 2*size(F,1));
ishift = -round((4.65 - ppmcen) * 2*size(F,1) * dt * B0);
for m = 1:size(F,2)
    x = [F(:, m); zeros(size(F,1), 1)];
    if conv == 3, d = fft(conj(x)); else, d = fft(x); end
    fprintf(fid, ' $NMUSED\n XTRASH = 0\n $END\n');
    fprintf(fid, ' $BASIS\n ID = ''%s'',\n METABO = ''%s'',\n CONC = 1.,\n TRAMP = 1.,\n VOLUME = 1.,\n ISHIFT = %d\n $END\n', ...
            names{m}, names{m}, ishift);
    fprintf(fid, '%15.6E%15.6E\n', [real(d) imag(d)].');
end
fclose(fid);
end

function write_raw(fname, x, conjflag, echot, B0)
if conjflag, x = conj(x); end
fid = fopen(fname, 'w');
fprintf(fid, ' $SEQPAR\n ECHOT = %.2f\n HZPPPM = %.5f\n SEQ = ''PULSEACQ''\n $END\n', echot, B0);
fprintf(fid, ' $NMID\n ID = ''real'', FMTDAT = ''(2E15.6)''\n VOLUME = 1., TRAMP = 1.\n $END\n');
fprintf(fid, '%15.6E%15.6E\n', [real(x(:)) imag(x(:))].');
fclose(fid);
end

function [tab, sd] = run_lcm(stem, fraw, fbas, dt, B0, N, ppmcen, echot, opts, chuse1, chcomb, key, owner, bin, ldir)
fc = [stem '.CONTROL']; ft = [stem '.table'];
fid = fopen(fc, 'w');
fprintf(fid, ' $LCMODL\n');
if ~isempty(key),   fprintf(fid, ' KEY = %d\n', key); end
if ~isempty(owner), fprintf(fid, ' OWNER = ''%s''\n', owner); end
fprintf(fid, ' TITLE = ''real 31P''\n');
fprintf(fid, ' FILRAW = ''%s''\n FILBAS = ''%s''\n', fraw, fbas);
fprintf(fid, ' FILTAB = ''%s''\n FILPS = ''%s.ps''\n FILCOO = ''%s.coord''\n', ft, stem, stem);
fprintf(fid, ' HZPPPM = %.5f\n DELTAT = %.6e\n NUNFIL = %d\n ECHOT = %.2f\n PPMCEN = %.4f\n', ...
        B0, dt, N, echot, ppmcen);
for k = 1:2:numel(opts)
    v = opts{k+1};
    if ischar(v)
        if any(strcmp(v, {'T','F'})), fprintf(fid, ' %s = %s\n', opts{k}, v);
        else, fprintf(fid, ' %s = ''%s''\n', opts{k}, v); end
    else
        fprintf(fid, ' %s = %g\n', opts{k}, v);
    end
end
fprintf(fid, ' NUSE1 = %d\n', numel(chuse1));
for k = 1:numel(chuse1), fprintf(fid, ' CHUSE1(%d) = ''%s''\n', k, chuse1{k}); end
fprintf(fid, ' NCOMBI = %d\n', numel(chcomb));
for k = 1:numel(chcomb), fprintf(fid, ' CHCOMB(%d) = ''%s''\n', k, chcomb{k}); end
fprintf(fid, ' $END\n');
fclose(fid);
here = pwd;
if ~isempty(ldir), cd(ldir); end
[st, out] = system(sprintf('"%s" < "%s"', bin, fc));
cd(here);
if st ~= 0, warning('LCModel returned %d: %s', st, out); end
[tab, sd] = parse_table(ft);
end

function [T, SD] = parse_table(ft)
T = containers.Map(); SD = containers.Map();
if ~exist(ft, 'file'), return; end
L = splitlines(fileread(ft));
i0 = find(startsWith(strtrim(L), '$$CONC'), 1);
if isempty(i0), return; end
for i = i0+2:numel(L)
    tok = regexp(L{i}, '^\s*([-+0-9.Ee]+)\s+(\d+)%\s+([-+0-9.Ee]+)\s+(\S+)\s*$', 'tokens', 'once');
    if isempty(tok)
        if startsWith(strtrim(L{i}), '$$'), break; end
        continue;
    end
    T(tok{4}) = str2double(tok{1}); SD(tok{4}) = str2double(tok{2});
end
end

function C = parse_coord(fc)
% Curves from an LCModel .coord file (LCOORD = 9).
C = struct('ppm', [], 'data', [], 'fit', [], 'bg', [], 'mets', {{}}, 'curves', []);
if ~exist(fc, 'file'), warning('No coord file %s', fc); return; end
L = splitlines(fileread(fc));
i = 1; NY = 0;
while i <= numel(L)
    s = L{i};
    if contains(s, 'points on ppm-axis')
        NY = str2double(regexp(s, '\d+', 'match', 'once')); [C.ppm, i] = readn(L, i+1, NY);
    elseif NY > 0 && contains(s, 'phased data points follow'), [C.data, i] = readn(L, i+1, NY);
    elseif NY > 0 && contains(s, 'points of the fit to the data follow'), [C.fit, i] = readn(L, i+1, NY);
    elseif NY > 0 && contains(s, 'background values follow'), [C.bg, i] = readn(L, i+1, NY);
    elseif NY > 0 && contains(s, 'Conc. =')
        nm = strtok(strtrim(s)); [v, i] = readn(L, i+1, NY);
        C.mets{end+1} = nm; C.curves(:, end+1) = v(:);
    else
        i = i + 1;
    end
end
if isempty(C.bg), C.bg = zeros(size(C.ppm)); end
if isempty(C.curves), C.curves = zeros(numel(C.ppm), 0); end
end

function [v, i] = readn(L, i, n)
v = zeros(n, 1); k = 0;
while k < n && i <= numel(L)
    x = sscanf(regexprep(L{i}, '[dD]([-+]?\d)', 'E$1'), '%f');
    x = x(1:min(numel(x), n - k));
    v(k+1:k+numel(x)) = x; k = k + numel(x); i = i + 1;
end
end

function v = report_amounts(amp, b, mets)
REP = {'PCr','Pi_intra','PE','PC','GPC','GPE','gATP','aATP','NADplus','NADH','tNAD','UDP'};
nP = containers.Map({'NADH','NADplus','UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}, {2,2,2,2,2,2,2});
c = @(m) conc1(amp, b, mets, m, nP);
v = nan(1, numel(REP));
for j = 1:numel(REP)
    switch REP{j}
        case 'tNAD', v(j) = c('NADH') + c('NADplus');
        case 'UDP'
            u = 0;
            for m = {'UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}
                if any(strcmp(mets, m{1})), u = u + c(m{1}); end
            end
            v(j) = u;
        otherwise, v(j) = c(REP{j});
    end
end
end

function x = conc1(amp, b, mets, m, nP)
i = strcmp(mets, m);
if ~any(i), x = NaN; return; end
S = sum(b.prior.rel_amp(strcmp(b.prior.metabolite, m)));
if isKey(nP, m), n = nP(m); else, n = 1; end
x = amp(i) * S / n;
end

function W = rep_matrix(b, mets, REP)
nP = containers.Map({'NADH','NADplus','UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}, {2,2,2,2,2,2,2});
udp = {'UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
W = zeros(numel(REP), numel(mets));
for m = 1:numel(mets)
    Ssum = sum(b.prior.rel_amp(strcmp(b.prior.metabolite, mets{m})));
    if isKey(nP, mets{m}), fac = Ssum / nP(mets{m}); else, fac = Ssum; end
    for j = 1:numel(REP)
        switch REP{j}
            case 'tNAD', hit = any(strcmp(mets{m}, {'NADH','NADplus'}));
            case 'UDP',  hit = any(strcmp(mets{m}, udp));
            otherwise,   hit = strcmp(mets{m}, REP{j});
        end
        if hit, W(j, m) = fac; end
    end
end
for j = 1:numel(REP)
    if ~any(W(j,:)) && ~any(strcmp(REP{j}, {'UDP','tNAD'})), W(j,:) = NaN; end
end
end

function p = wsl_path(p)
if numel(p) > 2 && p(2) == ':'
    p = ['/mnt/' lower(p(1)) strrep(p(3:end), '\', '/')];
end
end
