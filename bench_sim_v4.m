%% BENCH_SIM_V4  Simulated spectra for Figs. 2-3 (no fitting here).
% Same truth as bench_nad_sim_v3 (amplitudes, in vivo shifts and widths from
% the fit of the real spectrum), two cases:
%   1 'no baseline'   : model lineshapes with the in vivo shifts and widths,
%                       nominal phase, flat baseline
%   2 '+ baseline'    : the same + broad hump (area 3 x PCr, 500 Hz at -2 ppm)
%                       and the first 2 FID points zeroed
% x NAD+/NADH ratios (RATIOS) x N_REP noise realisations x NOISE_FACTORS
% (x the real single-voxel noise SD). The noise realisation of repetition k
% is the same at every level and in both cases (paired), so only the amount
% of noise or the baseline changes. NOISE_FACTORS(1) = 0.01 is "noiseless"
% (PCr SNR ~10^4); exactly zero noise is avoided because the fitters estimate
% the noise level from the data.
%
% Output (same layout as bench_nad_sim_v3.mat, EST empty), so the fitter
% scripts run on it unchanged:
%   bench_refit       BENCH_FILE = 'bench_sim_v4.mat'; OUTFILE = 'refit_v4.mat';
%                     VARIANTS = {'NEW K20 reg1 eta',20,1,true,'NEW'; 'OLD K20 reg1 eta',20,1,true,'OLD'};
%   lcm_benchmark     BENCH_FILE = 'bench_sim_v4.mat'; OUTDIR = fullfile(pwd,'lcm_v4'); RESULTS = 'lcm_v4.mat';
%   amares_benchmark  BENCH_FILE = 'bench_sim_v4.mat'; RESULTS = 'amares_v4.mat';
%   fsl_export        BENCH_FILE = 'bench_sim_v4.mat'; OUTDIR = fullfile(pwd,'fsl_v4'); then in WSL
%                     python fsl_benchmark.py .../fsl_v4/fsl_input.mat --scen 1 2 \
%                         --configs "FSL-NEW spline grpA" "FSL-OLD spline grpA" --workers 8
% Needs: basis_real in the workspace, make_basis_ren2020, bench_nad_sim_v3.mat (truth).

SRC_FILE      = 'bench_nad_sim_v3.mat';
OUTFILE       = 'bench_sim_v4.mat';
NOISE_FACTORS = [0.01 0.25 0.5 1 2 3 5];   % PCr SNR ~ 10700, 430, 215, 107, 54, 36, 21
N_REP         = 12;
HUMP_AREA     = 3;  HUMP_PPM = -2;  HUMP_LW = 500;  ROLL_PTS = 2;
SEED0         = 4000;

S0 = load(SRC_FILE, 'ampT', 'nuis', 'sigma0', 'T0_TRUE', 'RATIOS', 'REP');
bTrue = make_basis_ren2020(basis_real, 'normalise_all', true);
B0 = bTrue.B0_MHz; ppm = bTrue.ppm(:); N = numel(ppm); metsT = bTrue.metabolites;
iH = strcmp(metsT, 'NADH'); iP = strcmp(metsT, 'NADplus'); iC = strcmp(metsT, 'PCr');
ampT = S0.ampT(:); tNAD = ampT(iH) + ampT(iP); RATIOS = S0.RATIOS; REP = S0.REP;
sigma0 = S0.sigma0; T0_TRUE = S0.T0_TRUE;

SCN = {'no baseline', '+ baseline'};
nS = numel(SCN); nQ = numel(RATIOS); nN = numel(NOISE_FACTORS); nR = numel(REP);
SPEC = complex(nan(N, nS, nQ, N_REP, nN)); TRU = nan(nQ, nR);
CLEAN = complex(nan(N, nS, nQ));
for q = 1:nQ
    a = ampT; a(iH) = tNAD/(1+RATIOS(q)); a(iP) = tNAD*RATIOS(q)/(1+RATIOS(q));
    TRU(q,:) = report_amounts(a, bTrue, metsT);
    for s = 1:nS
        o = struct('t0', T0_TRUE, 'hump', HUMP_AREA*(s == 2), 'roll', ROLL_PTS*(s == 2));
        CLEAN(:, s, q) = simulate(bTrue, a, S0.nuis, o, HUMP_PPM, HUMP_LW);
        for k = 1:N_REP
            rng(SEED0 + k);
            nz = randn(N, 1) + 1i*randn(N, 1);
            for n = 1:nN
                SPEC(:, s, q, k, n) = CLEAN(:, s, q) + sigma0 * NOISE_FACTORS(n) * nz;
            end
        end
    end
end
pk = max(real(CLEAN(abs(ppm) < 0.5, 1, 1) .* exp(-1i*2*pi*B0*T0_TRUE*ppm(abs(ppm) < 0.5))));
SNR_PCR = pk / sigma0 ./ NOISE_FACTORS;
fprintf('PCr SNR per noise level: %s\n', mat2str(round(SNR_PCR)));

out = struct('SPEC', SPEC, 'TRU', TRU, 'REP', {REP}, 'SCN', {SCN}, 'RATIOS', RATIOS, ...
             'NOISE_FACTORS', NOISE_FACTORS, 'ampT', ampT, 'nuis', S0.nuis, 'sigma0', sigma0, ...
             'T0_TRUE', T0_TRUE, 'ppm', ppm, 'SNR_PCR', SNR_PCR, 'CLEAN', CLEAN);
out.cfgn = {};
out.EST  = nan(nS, nQ, N_REP, 0, nR, nN);
save(OUTFILE, '-struct', 'out', '-v7.3');
fprintf('Saved %s: %d cases x %d ratios x %d reps x %d noise levels = %d spectra\n', ...
    OUTFILE, nS, nQ, N_REP, nN, nS*nQ*N_REP*nN);

figure('Color','w', 'Name', 'bench_sim_v4: examples');
show = [1 4 nN];
for i = 1:numel(show)
    subplot(numel(show), 1, i); hold on;
    for s = 1:nS
        y = real(SPEC(:, s, 1, 1, show(i)) .* exp(-1i*2*pi*B0*T0_TRUE*ppm));
        plot(ppm, y + (s-1)*0.3*pk, 'LineWidth', 0.7);
    end
    set(gca, 'XDir','reverse'); xlim([-20 10]); box off;
    title(sprintf('noise x%g (PCr SNR %.0f): %s (offset)', NOISE_FACTORS(show(i)), SNR_PCR(show(i)), strjoin(SCN, ' / ')));
end

%% ===========================================================================
function spec = simulate(b, amp, nuis, o, hump_ppm, hump_lw)
% As in bench_nad_sim_v3 (scenario B lineshapes: in vivo group shifts and widths).
t = b.t(:); B0 = b.B0_MHz; mets = b.metabolites;
gmap = {{'PE','PC'}, {'GPC','GPE'}, {'Pi_intra','Pi_extra'}, {'NADH','NADplus'}, {'aATP'}, {'gATP'}};
fid = zeros(numel(t), 1);
for m = 1:numel(mets)
    if amp(m) == 0, continue; end
    fsh = 0; lwo = 0;
    for g = 1:numel(gmap)
        if any(strcmp(gmap{g}, mets{m})), fsh = fsh + nuis.grp(g); end
    end
    jm = strcmp(nuis.mets, mets{m});
    if any(jm), lwo = nuis.lwoff(jm); end
    pk = b.prior(strcmp(b.prior.metabolite, mets{m}), :);
    for r = 1:height(pk)
        lw = max(pk.lw_T2_hz(r) + nuis.lw_inh + lwo, 0.5);
        fid = fid + amp(m) * comp(t, (pk.cs_ppm(r) + fsh)*B0, lw, pk.rel_amp(r), b.eta);
    end
end
if o.hump > 0
    fid = fid + o.hump * amp(strcmp(mets, 'PCr')) * comp(t, hump_ppm*B0, hump_lw, 1, 0);
end
if o.roll > 0, fid(1:o.roll) = 0; end
spec = fftshift(fft(fid)) .* exp(1i*2*pi*B0*o.t0*b.ppm(:));
end

function s = comp(t, cs_Hz, lw_hz, a, eta)
R2s = pi*lw_hz; sg = sqrt(2*log(2)) / (pi*lw_hz + eps);
s = a .* exp(1i*2*pi*cs_Hz.*t) .* ((1-eta)*exp(-R2s.*t) + eta*exp(-t.^2./(2*sg^2)));
end

function v = report_amounts(amp, b, mets)
REP = {'PCr','Pi_intra','PE','PC','GPC','GPE','gATP','aATP','NADplus','NADH','tNAD','UDP'};
nP = containers.Map({'NADH','NADplus','UDPG'}, {2,2,2});
v = nan(1, numel(REP));
for j = 1:numel(REP)
    switch REP{j}
        case 'tNAD', v(j) = conc1(amp, b, mets, 'NADH', nP) + conc1(amp, b, mets, 'NADplus', nP);
        case 'UDP',  v(j) = conc1(amp, b, mets, 'UDPG', nP);
        otherwise,   v(j) = conc1(amp, b, mets, REP{j}, nP);
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
