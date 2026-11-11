%% FSL_EXPORT
% Export the simulated benchmark spectra and both basis sets for the FSL-MRS
% arm of the benchmark (fsl_benchmark.py, run in WSL). Same inputs as the
% LCModel arm: same spectra (nominal t0 phase removed, scaled so PCr = 1),
% same basis FIDs (Lorentzian lw_T2 per line, MP_broad dropped).
%
% Writes OUTDIR/fsl_input.mat (-v7, readable by scipy). Then in WSL:
%   python fsl_benchmark.py /mnt/c/.../fsl_bench/fsl_input.mat --workers 8
% and add fullfile('fsl_bench','fsl_results.mat') to LCM_FILES in bench_compare.
%
% Needs in the workspace: basis_real. Files: bench_nad_sim_v3.mat, make_basis_ren2020.m.

BENCH_FILE = 'bench_sim_v4.mat';
OUTDIR     = fullfile(pwd, 'fsl_bench_v4');
LW_CAL     = 8;      % Hz extra Lorentzian broadening of the calibration spectrum (as for LCModel)

if ~exist(OUTDIR, 'dir'), mkdir(OUTDIR); end
S = load(BENCH_FILE);
bOld = basis_real;
bNew = make_basis_ren2020(basis_real, 'normalise_all', true);
B0 = bNew.B0_MHz; t = bNew.t(:); N = numel(t); dt = t(2) - t(1); ppm = bNew.ppm(:);

% axis check: FSL-MRS uses ppm = f/B0 for 31P (no shift), so the carrier must be 0 ppm
f = ((0:N-1)' - floor(N/2)) / (N*dt);
pfit = polyfit(f/B0, ppm, 1);
fprintf('ppm axis: slope %.4f (expect +1), carrier %.3f ppm (expect 0)\n', pfit(1), pfit(2));
if abs(pfit(1) - 1) > 1e-3, error('ppm axis is not +f/B0 + const.'); end
if abs(pfit(2)) > 0.05, warning('Carrier is not 0 ppm: FSL-MRS ppm axis (and ppmlim) will be offset by %.3f ppm.', pfit(2)); end

DSCALE = 1 / S.ampT(strcmp(bNew.metabolites, 'PCr'));

% spectra -> FIDs, exactly as in lcm_benchmark (before its RAW conjugation)
[~, nS, nQ, nK, nN] = size(S.SPEC);
phi1nom = 2*pi*B0*S.T0_TRUE;
FID = complex(zeros(N, nS, nQ, nK, nN));
for s = 1:nS, for q = 1:nQ, for k = 1:nK, for n = 1:nN
    spec = S.SPEC(:, s, q, k, n) .* exp(-1i*phi1nom*ppm);
    FID(:, s, q, k, n) = ifft(ifftshift(spec)) * DSCALE;
end, end, end, end

% basis FIDs, REP mapping matrices and calibration amplitudes
out = struct('FID', FID, 'dt', dt, 'B0', B0, 'DSCALE', DSCALE, 'LW_CAL', LW_CAL, ...
             'dims', [N nS nQ nK nN]);
out.REP = S.REP;
bs = {'OLD', bOld; 'NEW', bNew};
for i = 1:2
    b = bs{i,2};
    [F, mets] = basis_fids(b);
    keep = ~strcmp(mets, 'MP_broad');
    F = F(:, keep); mets = mets(keep);
    aCal = zeros(numel(mets), 1);
    for m = 1:numel(mets)
        j = strcmp(bNew.metabolites, mets{m});
        if any(j), aCal(m) = S.ampT(j) * DSCALE; end
    end
    out.(['F_' bs{i,1}])    = F;
    out.(['mets_' bs{i,1}]) = mets;
    out.(['W_' bs{i,1}])    = rep_matrix(b, mets, S.REP);
    out.(['aCal_' bs{i,1}]) = aCal;
    fprintf('%s basis: %d metabolites: %s\n', bs{i,1}, numel(mets), strjoin(mets, ' '));
end
fout = fullfile(OUTDIR, 'fsl_input.mat');
save(fout, '-struct', 'out', '-v7');
fprintf('Saved %s (%d spectra). WSL path: %s\n', fout, nS*nQ*nK*nN, wsl_path(fout));

%% ===========================================================================
function [F, mets] = basis_fids(b)
% Same as lcm_benchmark: Lorentzian lw_T2 per line, the fitter adds the rest.
t = b.t(:); B0 = b.B0_MHz; mets = b.metabolites(:)';
F = zeros(numel(t), numel(mets));
for m = 1:numel(mets)
    pk = b.prior(strcmp(b.prior.metabolite, mets{m}), :);
    for r = 1:height(pk)
        F(:, m) = F(:, m) + pk.rel_amp(r) * exp(1i*2*pi*pk.cs_ppm(r)*B0*t - pi*pk.lw_T2_hz(r)*t);
    end
end
end

function W = rep_matrix(b, mets, REP)
% EST = W * amp reproduces report_amounts (per molecule, 2 31P for NAD/UDP).
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
for j = 1:numel(REP)      % metabolite absent from the basis -> NaN (as report_amounts)
    if ~any(W(j,:)) && ~any(strcmp(REP{j}, {'UDP','tNAD'})), W(j,:) = NaN; end
end
end

function p = wsl_path(p)
if numel(p) > 2 && p(2) == ':'
    p = ['/mnt/' lower(p(1)) strrep(p(3:end), '\', '/')];
end
end
