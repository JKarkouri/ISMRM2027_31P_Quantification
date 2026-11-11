%% BENCH_REFIT
% Refit the SAME simulated spectra (bench_nad_sim_v3.mat) with variants of the
% in-house full model: number of baseline knots, regularised baseline
% (curvature prior) and fitted lineshape (eta). Results go to refit_results.mat
% and are merged with the earlier configurations and LCModel by bench_compare.
%
% Needs: basis_real, x0_real in the workspace; make_basis_ren2020 and the
% latest fit_31P_spectrum_v4_cs (options 'baseline_curv_sd', 'fit_eta').

BENCH_FILE = 'bench_sim_v4.mat';
OUTFILE    = 'refit_results_v4.mat';
SCEN_RUN   = [];                      % [] = all scenarios present in BENCH_FILE
% name, n_baseline_knots, baseline_curv_sd ([] = off), fit_eta, basis ('NEW' or 'OLD'; optional, default NEW)
VARIANTS = {'NEW K20 reg1 eta', 20, 1, true, 'NEW';
            'OLD K20 reg1 eta', 20, 1, true, 'OLD'};

S = load(BENCH_FILE);
bNew = make_basis_ren2020(basis_real, 'normalise_all', true);
LGW  = {'gATP', {'gATP'}, 15; 'aATP', {'aATP'}, 15; 'Pi', {'Pi_intra','Pi_extra'}, 15};
BB = struct('NEW', bNew, 'OLD', basis_real);
NADP = struct('NEW', [0 5], 'OLD', [15 10]);       % NAD offset prior vs alpha-ATP (as in bench_nad_sim)
if size(VARIANTS, 2) < 5, VARIANTS(:, 5) = {'NEW'}; end
BASE = struct();
for bn = {'NEW', 'OLD'}
    b = BB.(bn{1}); M = numel(b.metabolites);
    AB = [-0.1*ones(1,M); Inf(1,M)]; AB(:, strcmp(b.metabolites,'MP_broad')) = [0; 1e-12];
    nn = b.metabolites(ismember(b.metabolites, {'NADH','NADplus','UDPG'}));
    P = NaN(4,2); P(4,:) = NADP.(bn{1});
    BASE.(bn{1}) = {'ppm_range',[-20 10], 'x0',x0_real, 'plot',false, 't0',S.T0_TRUE, 't0_prior_sd',25e-6, ...
        'atp_shift_prior_hz',[0 8; 0 8], 'nad_rel_aatp',true, 'lw_groups', LGW, ...
        'amp_bounds', AB, 'nonneg_mets', nn, 'fsh_grp_prior_hz', P};
end

[~, nS, nQ, nK, nN] = size(S.SPEC);
if isempty(SCEN_RUN)
    has = any(any(any(~isnan(S.SPEC(1,:,:,:,:)), 3), 4), 5);
    SCEN_RUN = find(has(:))';
end
nV = size(VARIANTS,1); nR = numel(S.REP);
EST_R = nan(nS, nQ, nK, nV, nR, nN); ETA = nan(nS, nQ, nK, nV, nN);
cfgR = VARIANTS(:,1)';
tic; nfit = 0; ntot = numel(SCEN_RUN)*nQ*nK*nN*nV;
for s = SCEN_RUN
    for q = 1:nQ
        for k = 1:nK
            for n = 1:nN
                spec = S.SPEC(:, s, q, k, n);
                for v = 1:nV
                    bv = BB.(VARIANTS{v,5});
                    args = [BASE.(VARIANTS{v,5}), {'n_baseline_knots', VARIANTS{v,2}, 'fit_eta', VARIANTS{v,4}}];
                    if ~isempty(VARIANTS{v,3}), args = [args, {'baseline_curv_sd', VARIANTS{v,3}}]; end %#ok<AGROW>
                    try
                        [~, r] = evalc('fit_31P_spectrum_v4_cs(spec, bv, args{:})');
                        EST_R(s,q,k,v,:,n) = report_amounts(r.amplitudes_raw(:), bv, r.metabolites);
                        ETA(s,q,k,v,n) = r.eta;
                    catch ME
                        warning('%s, s%d q%d k%d: %s', cfgR{v}, s, q, k, ME.message);
                    end
                    nfit = nfit + 1;
                end
            end
        end
        fprintf('  %d/%d fits, %.0f s\n', nfit, ntot, toc);
    end
end
save(OUTFILE, 'EST_R', 'cfgR', 'ETA', 'VARIANTS', 'SCEN_RUN');
fprintf('Saved %s. Fitted eta (truth: %.2f in E/F, basis eta elsewhere): %s\n', OUTFILE, 0.3, ...
    mat2str(squeeze(mean(mean(mean(mean(ETA, 2, 'omitnan'), 3, 'omitnan'), 5, 'omitnan'), 1, 'omitnan'))', 3));

%% ===========================================================================
function v = report_amounts(amp, b, mets)
REP = {'PCr','Pi_intra','PE','PC','GPC','GPE','gATP','aATP','NADplus','NADH','tNAD','UDP'};
nP = containers.Map({'NADH','NADplus','UDPG'}, {2,2,2});
c = @(m) conc1(amp, b, mets, m, nP);
v = nan(1, numel(REP));
for j = 1:numel(REP)
    switch REP{j}
        case 'tNAD', v(j) = c('NADH') + c('NADplus');
        case 'UDP',  v(j) = c('UDPG');
        otherwise,   v(j) = c(REP{j});
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
