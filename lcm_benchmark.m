%% LCM_BENCHMARK
% Run LCModel on the SAME simulated spectra as bench_nad_sim (v3 file) with
% basis sets built from the SAME basis structs (old basis_real and the Ren
% basis), so LCModel and the in-house fitter can be compared in the 2x2
% design basis x fitter. Results are printed in the same format and appended
% as configurations 'LCM-OLD' and 'LCM-NEW'.
%
% Steps
%   1. Check the ppm/frequency axis of the basis struct.
%   2. Write one .BASIS per basis struct (LCModel names <= 6 characters).
%   3. CONVENTION CHECK: a noiseless spectrum built from the basis FIDs
%      themselves is fitted with every combination of the two unknown file
%      conventions (basis spectrum ordering/sign, RAW conjugation). The
%      combination that returns the known amplitudes is used for everything
%      (and the absolute LCModel scale is calibrated on it). If none works,
%      the script stops - nothing silently wrong.
%   4. Export every simulated spectrum as .RAW (nominal t0 phase removed, as
%      the in-house fit is also given the nominal t0), write .CONTROL with
%      CIBM-like settings, run LCModel, parse the .TABLE.
%   5. Bias tables vs truth, merged with the in-house configurations.
%
% Needs in the workspace: basis_real. Files: bench_nad_sim_v3.mat,
% make_basis_ren2020.m. LCModel binary path and licence below.

%% ---- settings ------------------------------------------------------------
BENCH_FILE = 'bench_sim_v4.mat';
OUTDIR     = fullfile(pwd, 'lcm_v4');
LCM_BIN    = 'LCModel_exe\LCModel.exe';
LCM_DIR    = '';
KEY        = [];
OWNER      = '';
ECHOT      = 0.35;
BASIS_LW   = 'T2';
SCEN_RUN   = [];                          % both cases of bench_sim_v4 (no baseline, + baseline)
RUNS       = {'OLD', 2; 'NEW', 2};        % both bases, 2 ppm knot spacing
RESULTS    = 'lcm_v4.mat';

% CIBM-like LCModel settings (Kaiser et al. 2025 supplement), adapted to this window
LCM_OPTS = { ...
    'PPMST', 10, 'PPMEND', -20, 'NSIMUL', 0, 'NRATIO', 0, 'DKNTMN', 2, 'RFWHM', 3, ...
    'SDDEGZ', 2, 'SDDEGP', 0.5, 'DEGZER', 0, 'DEGPPM', 0, 'DESDT2', 10, 'DEEXT2', 10, ...
    'DESDSH', 0.01, 'DOECC', 'F', 'DOWS', 'F', 'DOREFS', 'F', 'NCALIB', 0, 'NEACH', 999, ...
    'LCOORD', 9, 'LTABLE', 7, 'NAMREL', 'PCr', 'CONREL', 1};
CHUSE1 = {'PCr','PiIn','aATP','gATP'};      % primary analysis on the big peaks
CHCOMB = {'NADH+NADp'};                      % also report tNAD (parsed separately)

if isempty(LCM_BIN), error('Set LCM_BIN to your LCModel executable.'); end
if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end

%% ---- load benchmark + bases ----------------------------------------------
S = load(BENCH_FILE);
bOld = basis_real;
bNew = make_basis_ren2020(basis_real, 'normalise_all', true);
B0 = bNew.B0_MHz; t = bNew.t(:); N = numel(t); dt = t(2) - t(1);
ppm = bNew.ppm(:);

% 1. axis check: ppm must be linear in the FFT frequency index
f = ((0:N-1)' - floor(N/2)) / (N*dt);            % fftshift frequency axis (Hz)
pfit = polyfit(f/B0, ppm, 1);
fprintf('ppm axis: slope %.4f (expect +1), carrier at %.3f ppm, max dev %.2g ppm\n', ...
    pfit(1), pfit(2), max(abs(polyval(pfit, f/B0) - ppm)));
if abs(pfit(1) - 1) > 1e-3
    error('ppm axis is not +f/B0 + const; check the basis ppm convention before exporting.');
end
PPMCEN = pfit(2);

% name maps (LCModel names kept short)
MAP = {'PCr','PCr'; 'Pi_intra','PiIn'; 'Pi_extra','PiEx'; 'PE','PE'; 'PC','PC'; ...
       'GPC','GPC'; 'GPE','GPE'; 'gATP','gATP'; 'aATP','aATP'; 'bATP','bATP'; ...
       'NADH','NADH'; 'NADplus','NADp'; 'UDPG','UDPG'; 'MP_broad','MPbr'; ...
       'UDPGal','UGal'; 'UDPGlc','UGlc'; 'UDPGalNAc','UGaN'; 'UDPGlcNAc','UGlN'};
short = @(m) MAP{strcmp(MAP(:,1), m), 2};

%% 2. basis FIDs (time domain) for each basis struct
BAS = struct('name', {'OLD','NEW'}, 'b', {bOld, bNew});
for i = 1:2
    [BAS(i).fids, BAS(i).mets] = basis_fids(BAS(i).b, BASIS_LW);
    keep = ~strcmp(BAS(i).mets, 'MP_broad');      % MP_broad is switched off in all fits
    BAS(i).fids = BAS(i).fids(:, keep); BAS(i).mets = BAS(i).mets(keep);
    BAS(i).short = cellfun(short, BAS(i).mets, 'UniformOutput', false);
end
fwhmba = mean(bNew.prior.lw_T2_hz) / B0;

%% 3. convention check on a noiseless spectrum made from the NEW basis FIDs
aT = S.ampT(:); metsT = bNew.metabolites;
DSCALE = 1 / aT(strcmp(metsT,'PCr'));   % data scaled so PCr amplitude = 1 (LCModel works in single precision)
aV = zeros(numel(BAS(2).mets), 1);
for m = 1:numel(BAS(2).mets), aV(m) = aT(strcmp(metsT, BAS(2).mets{m})); end
fidV = BAS(2).fids * (aV*DSCALE) .* exp(-pi*8*t);  % + 8 Hz extra Lorentzian broadening
rng(1); fidV = fidV + 1e-3*max(abs(fidV))*(randn(N,1) + 1i*randn(N,1));   % tiny noise: LCModel needs a noise estimate
vdir = fullfile(OUTDIR, 'convention_check'); if ~exist(vdir, 'dir'), mkdir(vdir); end
best = struct('err', Inf);
fprintf('\nConvention check (near-noiseless, NEW basis):\n');
% Convention verified against LCModel 6.3-1N (open-source build): RAW = conj of
% our FID, basis = FFT of that (conv 3), with ISHIFT compensating LCModel's
% assumption that basis spectra are centred at 4.65 ppm. Conv 1 / no conj is
% the mirror-image equivalent, tried as a fallback.
COMBOS = [3 1; 1 0];
for ic = 1:size(COMBOS,1)
    bc = COMBOS(ic,1); rc = COMBOS(ic,2);
    if true
        fb = fullfile(vdir, sprintf('chk_b%d.BASIS', bc));
        write_basis(fb, BAS(2).fids, BAS(2).short, dt, B0, ECHOT, fwhmba, bc, PPMCEN);
        fr = fullfile(vdir, sprintf('chk_b%d_r%d.RAW', bc, rc));
        write_raw(fr, fidV, rc, ECHOT, B0);
        tab = run_lcm(fullfile(vdir, sprintf('chk_b%d_r%d', bc, rc)), fr, fb, dt, B0, N, PPMCEN, ...
                      ECHOT, LCM_OPTS, CHUSE1, {}, KEY, OWNER, LCM_BIN, LCM_DIR);
        c = conc_vector(tab, BAS(2).short);
        nz = sum(c > 0); iPc = strcmp(BAS(2).mets, 'PCr');
        big = ismember(BAS(2).mets, {'PCr','Pi_intra','gATP','aATP','GPC','PE','GPE'}) & aV' > 0;
        iP = strcmp(BAS(2).mets, 'PCr');
        rel_est = c(big) / c(iP); rel_tru = aV(big)' / aV(iP);
        err = sqrt(mean(log(max(rel_est, 1e-6) ./ rel_tru).^2));
        if ~isfinite(err), err = Inf; end
        fprintf('  basis convention %d, RAW conj %d: rms log-ratio error %.3f | PCr conc %.3g | %d/%d metabolites > 0 | %s\n', ...
            bc, rc, err, c(iPc), nz, numel(c), table_msg(fullfile(vdir, sprintf('chk_b%d_r%d.table', bc, rc))));
        if err < best.err, best = struct('err', err, 'bc', bc, 'rc', rc, 'scale', c(iP)/(aV(iP)*DSCALE)); end
    end
    if best.err < 0.05, break; end
end
if best.err > 0.05
    error(['No file convention reproduced the noiseless amplitudes (best rms error %.3f). ' ...
           'Open the .ps files in %s to see what LCModel did.'], best.err, vdir);
end
fprintf('Using basis convention %d, RAW conj %d (error %.4f); LCModel scale = %.4g x amplitude\n', ...
    best.bc, best.rc, best.err, best.scale);

for i = 1:2
    BAS(i).file = fullfile(OUTDIR, sprintf('%s.BASIS', BAS(i).name));
    write_basis(BAS(i).file, BAS(i).fids, BAS(i).short, dt, B0, ECHOT, fwhmba, best.bc, PPMCEN);
end

%% 4. run LCModel on every simulated spectrum
[~, nS, nQ, nK, nN] = size(S.SPEC);
if isempty(SCEN_RUN), SCEN_RUN = 1:nS; end
nR = numel(S.REP);
nRun = size(RUNS,1);
EST_L = nan(nS, nQ, nK, nRun, nR, nN);
cfgL = cell(1, nRun);
for ir = 1:nRun, cfgL{ir} = sprintf('LCM-%s dk%g', RUNS{ir,1}, RUNS{ir,2}); end
phi1nom = 2*pi*B0*S.T0_TRUE;
tic; nrun = 0; ntot = numel(SCEN_RUN)*nQ*nK*nN*nRun;
for s = SCEN_RUN
    for q = 1:nQ
        for k = 1:nK
            for n = 1:nN
                spec = S.SPEC(:, s, q, k, n) .* exp(-1i*phi1nom*ppm);   % remove nominal t0 phase
                fid  = ifft(ifftshift(spec)) * DSCALE;
                tag  = sprintf('s%d_q%d_k%02d_n%d', s, q, k, n);
                fr   = fullfile(OUTDIR, 'raw', [tag '.RAW']);
                if ~exist(fileparts(fr), 'dir'), mkdir(fileparts(fr)); end
                write_raw(fr, fid, best.rc, ECHOT, B0);
                for ir = 1:nRun
                    i = find(strcmp({BAS.name}, RUNS{ir,1}));
                    opts = LCM_OPTS; opts{find(strcmp(opts(1:2:end), 'DKNTMN'))*2} = RUNS{ir,2};
                    od = fullfile(OUTDIR, sprintf('%s_dk%g', RUNS{ir,1}, RUNS{ir,2}));
                    if ~exist(od, 'dir'), mkdir(od); end
                    tab = run_lcm(fullfile(od, tag), fr, BAS(i).file, dt, B0, N, PPMCEN, ECHOT, ...
                                  opts, CHUSE1, CHCOMB, KEY, OWNER, LCM_BIN, LCM_DIR);
                    c = conc_vector(tab, BAS(i).short) / best.scale / DSCALE; % -> amplitude units
                    amp = zeros(numel(BAS(i).b.metabolites), 1);
                    for m = 1:numel(BAS(i).mets)
                        amp(strcmp(BAS(i).b.metabolites, BAS(i).mets{m})) = c(m);
                    end
                    EST_L(s,q,k,ir,:,n) = report_amounts(amp, BAS(i).b, BAS(i).b.metabolites);
                    nrun = nrun + 1;
                end
            end
        end
        fprintf('  %d/%d LCModel runs, %.0f s\n', nrun, ntot, toc);
    end
end
save(fullfile(OUTDIR, RESULTS), 'EST_L', 'cfgL', 'best', 'SCEN_RUN', 'RUNS');

%% 5. tables: in-house configurations + LCModel
EST = cat(4, S.EST, EST_L);
cfgn = [S.cfgn, cfgL];
REP = S.REP; TRU = S.TRU; pcr_true = TRU(1, strcmp(REP,'PCr'));
n0 = find(S.NOISE_FACTORS == 1, 1); if isempty(n0), n0 = 1; end
for q = 1:nQ
    fprintf('\n===== true NAD+/NADH = %g, noise x%g =====\n', S.RATIOS(q), S.NOISE_FACTORS(n0));
    for s = SCEN_RUN
        fprintf('-- %s\n', S.SCN{s});
        print_nad_rows(EST, TRU, REP, cfgn, s, q, n0, pcr_true);
    end
end
for s = SCEN_RUN(SCEN_RUN >= 4)
    fprintf('\n===== All metabolites, %s, ratios pooled: bias %% / RMSE %% =====\n', S.SCN{s});
    fprintf('%-9s', ''); fprintf(' | %16s', cfgn{:}); fprintf('\n');
    for j = 1:nR
        fprintf('%-9s', REP{j});
        for c = 1:numel(cfgn)
            [b, rm] = pooled(EST, TRU, s, c, j, n0);
            fprintf(' | %+6.1f / %6.1f ', b, rm);
        end
        fprintf('\n');
    end
end

%% ===========================================================================
function [F, mets] = basis_fids(b, lwmode)
% One time-domain FID per metabolite, same line model as the fitter (Lorentzian
% part only, lw_T2 per line; LCModel adds its own broadening).
t = b.t(:); B0 = b.B0_MHz; mets = b.metabolites(:)';
F = zeros(numel(t), numel(mets));
for m = 1:numel(mets)
    pk = b.prior(strcmp(b.prior.metabolite, mets{m}), :);
    for r = 1:height(pk)
        lw = pk.lw_T2_hz(r); if ~strcmp(lwmode, 'T2'), lw = 1; end
        F(:, m) = F(:, m) + pk.rel_amp(r) * exp(1i*2*pi*pk.cs_ppm(r)*B0*t - pi*lw*t);
    end
end
end

function write_basis(fname, F, names, dt, B0, echot, fwhmba, conv, ppmcen)
% LCModel assumes basis spectra have their zero frequency at 4.65 ppm and
% shifts them by nint((4.65 - PPMCEN)/ppminc_basis) points; ISHIFT cancels that.
% .BASIS file. 'conv' selects how the frequency-domain basis is stored; the
% correct one is found by the convention check.
fid = fopen(fname, 'w');
fprintf(fid, ' $SEQPAR\n FWHMBA = %.5f,\n HZPPPM = %.5f,\n ECHOT = %.2f,\n SEQ = ''PULSEACQ''\n $END\n', ...
        fwhmba, B0, echot);
fprintf(fid, ' $BASIS1\n IDBASI = ''bench'',\n FMTBAS = ''(2E15.6)'',\n BADELT = %.6e,\n NDATAB = %d\n $END\n', ...
        dt, 2*size(F,1));
ishift = -round((4.65 - ppmcen) * 2*size(F,1) * dt * B0);
for m = 1:size(F,2)
    x = [F(:, m); zeros(size(F,1), 1)];      % LCModel basis: NDATAB = 2*NUNFIL (zero-filled)
    switch conv
        case 1, d = fft(x);
        case 2, d = fftshift(fft(x));
        case 3, d = fft(conj(x));
        case 4, d = fftshift(fft(conj(x)));
        case 5, d = flipud(fftshift(fft(x)));
        case 6, d = flipud(fftshift(fft(conj(x))));
    end
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
fprintf(fid, ' $NMID\n ID = ''sim'', FMTDAT = ''(2E15.6)''\n VOLUME = 1., TRAMP = 1.\n $END\n');
fprintf(fid, '%15.6E%15.6E\n', [real(x(:)) imag(x(:))].');
fclose(fid);
end

function m = table_msg(ft)
% First line of LCModel's error/diagnostic block, if any.
m = 'no .table';
if ~exist(ft, 'file'), return; end
txt = fileread(ft);
tok = regexp(txt, '(?m)^\s*(FATAL|ERROR|Warning|WARNING)[^\n]*', 'match', 'once');
if isempty(tok), m = 'table ok'; else, m = strtrim(tok); end
end

function tab = run_lcm(stem, fraw, fbas, dt, B0, N, ppmcen, echot, opts, chuse1, chcomb, key, owner, bin, ldir)
fc = [stem '.CONTROL']; ft = [stem '.table'];
fid = fopen(fc, 'w');
fprintf(fid, ' $LCMODL\n');
if ~isempty(key),   fprintf(fid, ' KEY = %d\n', key); end
if ~isempty(owner), fprintf(fid, ' OWNER = ''%s''\n', owner); end
fprintf(fid, ' TITLE = ''%s''\n', stem(max(1,end-30):end));
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
if st ~= 0, warning('LCModel returned %d for %s: %s', st, fc, out); end
fo = fopen([stem '.stdout.txt'], 'w'); fprintf(fo, '%s', out); fclose(fo);
tab = parse_table(ft);
end

function T = parse_table(ft)
% Returns a containers.Map name -> concentration from the $$CONC block.
T = containers.Map();
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
    T(tok{4}) = str2double(tok{1});
end
end

function c = conc_vector(T, names)
c = nan(1, numel(names));
for m = 1:numel(names)
    if isKey(T, names{m}), c(m) = T(names{m}); elseif ~isempty(T), c(m) = 0; end
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

function print_nad_rows(EST, TRU, REP, cfgn, s, q, n, pcr_true)
fprintf('%-10s | %16s %16s %16s %16s | %7s %7s\n', 'config', 'NAD+ bias%', 'NADH bias%', ...
    'tNAD bias%', 'UDP bias%', 'med R', 'NADH=0');
for c = 1:numel(cfgn)
    fprintf('%-10s |', cfgn{c});
    for nm = {'NADplus','NADH','tNAD','UDP'}
        j = strcmp(REP, nm{1});
        e = squeeze(EST(s,q,:,c,j,n)); e = e(~isnan(e)); tt = TRU(q,j);
        if isempty(e), fprintf(' %16s', 'n/a'); continue; end
        fprintf(' %+8.1f (%5.1f)', 100*(mean(e)-tt)/tt, 100*std(e)/sqrt(numel(e))/abs(tt));
    end
    eH = squeeze(EST(s,q,:,c,strcmp(REP,'NADH'),n)); eP = squeeze(EST(s,q,:,c,strcmp(REP,'NADplus'),n));
    ok = ~isnan(eH) & ~isnan(eP); eH = eH(ok); eP = eP(ok);
    z = eH < 1e-6*pcr_true; R = eP(~z) ./ eH(~z);
    if isempty(R), mr = Inf; else, mr = median(R); end
    fprintf(' | %7.1f %6.0f%%\n', mr, 100*mean(z));
end
end

function [b, rm] = pooled(EST, TRU, s, c, j, n)
e = []; tt = [];
for q = 1:size(TRU,1)
    ee = squeeze(EST(s,q,:,c,j,n)); e = [e; ee]; tt = [tt; TRU(q,j)*ones(size(ee))]; %#ok<AGROW>
end
ok = ~isnan(e);
b = 100*mean((e(ok) - tt(ok))./tt(ok)); rm = 100*sqrt(mean(((e(ok) - tt(ok))./tt(ok)).^2));
end
