function result = fit_31P_spectrum_v4_cs(spec_obs, basis, varargin)
% FIT_31P_SPECTRUM  Frequency-domain least-squares quantification of 31P MRS.
%
% Signal model:
%   Y(v) = B(v) + exp(i*(phi0 + v*phi1)) * sum_m C_m * M_m(v; fsh, lw_inh)
%
%   B(v) = natural cubic spline with n_knots uniformly-spaced knots
%          Columns are thin-plate radial basis functions r^3 + linear trend.
%          This is equivalent to a smoothing spline and far more stable than
%          high-order polynomials over wide windows (LCModel uses ~20 knots).
%
% Nonlinear parameters (3, optimised with lsqnonlin):
%   phi0, freq_shift, lw_inhom_hz
%
% Linear parameters (M + n_knots+2, solved analytically via lsqlin):
%   metabolite amplitudes C_1..C_M  +  spline baseline coefficients
%
% t0 is fixed (pass via 't0' parameter — not optimised).
%
% USAGE
%   result = fit_31P_spectrum(spec_obs, basis)
%   result = fit_31P_spectrum(spec_obs, basis, 'ppm_range',[-20 10], ...
%               'n_baseline_knots', 10, 'plot', true)
%
% KEY OPTIONAL PARAMETERS
%   'n_baseline_knots' : number of spline knots                 default 10
%                        4  = very smooth (similar to polynomial order 2)
%                        10 = recommended for wide 31P window
%                        20 = LCModel-like flexibility
%   'ppm_range'        : [lo hi] fitting window (ppm)           default [-20 10]
%   'plot'             : true/false                              default false
%   'lw_inhom_bounds'  : [lo hi] Hz                             default [0 100]
%   'use_ratio_constraints' : physiological amplitude bounds     default false
%   'ratio_constraints' : custom hard ratio bounds, N x 4 cell
%                     {num, den, lo, hi}: lo <= num/den <= hi (Inf/NaN hi = no
%                     upper bound). Replaces the built-in set when given,
%                     e.g. {'NADplus','NADH',5,Inf}.            default {}
%
% OPTIONS ADDED 2026-10-02 (defaults reproduce the previous behaviour exactly)
%   'fsh_grp_bound' : bound (ppm) on the per-group frequency offsets, order
%                     {PME, PDE, Pi, NAD}. Scalar (all groups) or 4-vector.
%                     0 for a group = offset fixed at 0 (removed from the
%                     nonlinear fit).                            default 0.3
%                     e.g. [0.3 0.3 0.3 0] fixes NAD to the global shift.
%   'clip_negative' : true  -> negative amplitudes set to 0 in
%                     result.amplitudes (previous behaviour)     default true
%                     result.amplitudes_raw is ALWAYS the unclipped solution.
%   'fit_imag'      : true  -> fit real AND imaginary parts jointly, with a
%                     separate spline baseline for the imaginary channel.
%                     result.baseline_imag_full holds that baseline. Metrics
%                     (CRLB etc.) are still computed on the real part.
%                                                                default false
%   'local_phase_mets' : cell array of metabolite names that share ONE extra
%                     zero-order phase offset (fitted), on top of phi0+phi1.
%                     Intended use: {'aATP','NADH','NADplus','UDPG'} to absorb
%                     a residual first-order phase error in the -7..-10 ppm
%                     region without fitting t0. Empty = off.   default {}
%   'local_phase_bound': bound (rad) on that offset              default 0.5
%
% OPTIONS ADDED 2026-10-07 (LCModel-style soft priors; defaults = off)
%   'nonneg_mets'   : cell of metabolite names whose amplitude has a hard
%                     lower bound of 0 INSIDE the fit (not clipping after),
%                     e.g. {'NADH','NADplus'}.                   default {}
%   'fsh_grp_prior_hz': [4 x 2] = [mean sd] in Hz for the per-group offsets
%                     {PME, PDE, Pi, NAD}; NaN row = no prior. Adds the
%                     penalty ((offset-mean)/sd)^2 (in noise units) so the
%                     offset can move only if the data clearly ask for it.
%                     e.g. [NaN NaN; NaN NaN; NaN NaN; 22 3]     default []
%   'atp_shift_prior_hz': [2 x 2] = [mean sd] in Hz for own shifts of
%                     {aATP; gATP} relative to the global (PCr) shift. A
%                     NaN row = that line keeps the global shift (default).
%                     Bound = |mean| + 4 SD. e.g. [0 8; NaN NaN] fits aATP.
%                     (Alternatively pass 6-element fsh_grp_bound / [6 x 2]
%                     fsh_grp_prior_hz in the order PME,PDE,Pi,NAD,aATP,gATP.)
%   'nad_rel_aatp'  : true -> the NAD group offset is measured from the aATP
%                     shift: NAD shift = global + aATP offset + NAD offset,
%                     and the NAD prior in fsh_grp_prior_hz applies to that
%                     RELATIVE offset. Requires an aATP shift.   default false
%                     result.fsh_groups always holds the EFFECTIVE offsets
%                     vs the global shift (so fsh_groups(4)*B0 is still the
%                     NAD offset vs the basis); result.nad_rel_aatp_hz is the
%                     fitted NAD-minus-aATP offset.
%   't0_prior_sd'   : SD (s) of a Gaussian prior on t0 centred on the 't0'
%                     value. If set, t0 (hence phi1) is FITTED, bounded to
%                     t0 +- 4 SD. e.g. 25e-6 for +-25 us.         default []
%
%   'phase_groups'  : AMARES-style regional phases. N x 3 cell, one row per
%                     group: {name, {metabolites}, prior_sd_deg}. Each group
%                     gets ONE extra zero-order phase (on top of phi0+phi1),
%                     with a Gaussian prior centred on 0 (bounded to +-4 SD).
%                     e.g. {'aATP_NAD', {'aATP','NADH','NADplus','UDPG'}, 10;
%                           'gATP', {'gATP'}, 10; 'PCr', {'PCr'}, 10}
%                                                                default {}
%   'lw_groups'     : regional linewidth offsets (Hz), same format:
%                     {name, {metabolites}, prior_sd_hz}; added to
%                     lw_inhom for those metabolites, prior centred on 0.
%                     e.g. {'gATP', {'gATP'}, 5; 'Pi', {'Pi_intra','Pi_extra'}, 5}
%                     Optional 4th column = prior MEAN (Hz for lw_groups, deg
%                     for phase_groups); bounds become mean +- 4 SD and the
%                     start value the mean, e.g. {'NADH', {'NADH'}, 5, 10}.
%                                                                default {}
%   'baseline_curv_sd': regularised baseline (LCModel-like). Gaussian prior on
%                     the second differences of the spline coefficients, SD
%                     in units of the noise SD. Smaller = stiffer. [] = off.
%                     e.g. 1 (stiff) ... 10 (mild).               default []
%   'fit_eta'       : true -> fit the Lorentz/Gauss mix eta (0..1) of all lines
%                     instead of using basis.eta.                 default false
%   'eta_prior'     : [mean sd] Gaussian prior on eta when fitted. default []
%   Fitted values: result.phase_group_deg, result.lw_group_hz, and per
%   metabolite result.met_phase (rad) / result.met_lw_offset (Hz).
%
% AIC/BIC now count ALL fitted parameters (amplitudes, baseline coefficients
% of every fitted channel, nonlinear parameters) and use the data residual
% only (prior penalties excluded). Values are therefore not comparable with
% AIC/BIC from earlier versions of this function.

%% ---- Parse --------------------------------------------------------------
p = inputParser;
addRequired(p,'spec_obs', @isnumeric);
addRequired(p,'basis',    @isstruct);
addParameter(p,'ppm_range',             [-20 10],  @(x)numel(x)==2);
addParameter(p,'n_baseline_knots',      10,        @(x)isscalar(x)&&x>=2&&x<=40);
addParameter(p,'plot',                  false,     @islogical);
addParameter(p,'x0',                    [],        @(x)isempty(x)||isstruct(x));
addParameter(p,'amp_bounds',            [],        @isnumeric);
addParameter(p,'lw_inhom_bounds',       [0 100],   @(x)numel(x)==2);
addParameter(p,'use_ratio_constraints', false,     @islogical);
addParameter(p,'ratio_constraints',     {},        @iscell);
addParameter(p,'t0',                    0,         @isnumeric);
addParameter(p,'fsh_grp_bound',         0.3,       @(x)isnumeric(x)&&all(x>=0));
addParameter(p,'clip_negative',         true,      @islogical);
addParameter(p,'fit_imag',              false,     @islogical);
addParameter(p,'local_phase_mets',      {},        @iscell);
addParameter(p,'local_phase_bound',     0.5,       @(x)isscalar(x)&&x>0);
addParameter(p,'nonneg_mets',           {},        @iscell);
addParameter(p,'fsh_grp_prior_hz',      [],        @isnumeric);
addParameter(p,'atp_shift_prior_hz',    [],        @isnumeric);
addParameter(p,'nad_rel_aatp',          false,     @islogical);
addParameter(p,'t0_prior_sd',           [],        @(x)isempty(x)||(isscalar(x)&&x>0));
addParameter(p,'phase_groups',          {},        @iscell);
addParameter(p,'lw_groups',             {},        @iscell);
addParameter(p,'baseline_curv_sd',      [],        @(x)isempty(x)||(isscalar(x)&&x>0));
addParameter(p,'fit_eta',               false,     @islogical);
addParameter(p,'eta_prior',             [],        @(x)isempty(x)||numel(x)==2);
parse(p, spec_obs, basis, varargin{:});

ppm_range        = sort(p.Results.ppm_range);
K                = p.Results.n_baseline_knots;
doPlot           = p.Results.plot;
x0_in            = p.Results.x0;
amp_bounds       = p.Results.amp_bounds;
lw_inhom_bounds  = p.Results.lw_inhom_bounds;
use_constraints  = p.Results.use_ratio_constraints;
ratio_cons_in    = p.Results.ratio_constraints;
t0_fixed         = p.Results.t0;
fsh_grp_bound_in = p.Results.fsh_grp_bound;
clip_negative    = p.Results.clip_negative;
fit_imag         = p.Results.fit_imag;
lph_mets         = p.Results.local_phase_mets;
lph_bound        = p.Results.local_phase_bound;
nonneg_mets      = p.Results.nonneg_mets;
grp_prior_in     = p.Results.fsh_grp_prior_hz;
atp_prior_in     = p.Results.atp_shift_prior_hz;
nad_rel_aatp     = p.Results.nad_rel_aatp;
t0_prior_sd      = p.Results.t0_prior_sd;
phase_groups_in  = p.Results.phase_groups;
lw_groups_in     = p.Results.lw_groups;
bl_curv_sd       = p.Results.baseline_curv_sd;
fit_eta          = p.Results.fit_eta;
eta_prior        = p.Results.eta_prior;

spec_obs = double(spec_obs(:));
ppm  = basis.ppm(:);
B0   = basis.B0_MHz;
eta  = basis.eta;
M    = numel(basis.metabolites);
N    = numel(ppm);
t    = basis.t(:);
fs   = basis.fs;

phi1_fixed = 2*pi * B0 * t0_fixed;
fprintf('  t0=%.4f ms -> phi1_fixed=%.4f rad/ppm\n', t0_fixed*1e3, phi1_fixed);

if numel(spec_obs) ~= N
    error('spec_obs has %d points but basis has %d.', numel(spec_obs), N);
end

%% ---- Fitting window -----------------------------------------------------
mask  = ppm >= ppm_range(1) & ppm <= ppm_range(2);
n_fit = sum(mask);
ppm_w = ppm(mask);

%% ---- Cubic B-spline baseline matrix -------------------------------------
% B-splines are LOCALLY SUPPORTED — each basis function is nonzero only
% near its knot. This is critical for wide 31P windows: regions with no
% signal (e.g. +10 to +20 ppm) contribute zero columns and don't distort
% the fit elsewhere. This is the same approach LCModel uses.
%
% Implementation: de Boor recursion on an augmented knot vector.
% n_bl = K+2 basis functions for K interior knots (order 4 = cubic).

% Augmented knot vector: clamp 4x at each boundary, K interior knots
t_int  = linspace(ppm_w(1), ppm_w(end), K+2);
t_aug  = [repmat(t_int(1),1,3), t_int, repmat(t_int(end),1,3)];
n_bl   = numel(t_aug) - 4;

P_w    = bspline_basis(ppm_w, t_aug, n_bl);
P_full = bspline_basis(ppm,   t_aug, n_bl);

%% ---- Noise estimate -----------------------------------------------------
noise_mask = ppm < -18 | ppm > 9;
if sum(noise_mask) > 20
    sigma_n = std(real(spec_obs(noise_mask)));
else
    r = real(spec_obs);
    sigma_n = median(abs(r - median(r))) / 0.6745;
end
sigma_n = max(sigma_n, 1e-12 * max(abs(real(spec_obs))));

%% ---- Normalise ----------------------------------------------------------
pcr_region = abs(ppm) < 2;
data_scale = max(abs(real(spec_obs(pcr_region))));
if data_scale < eps, data_scale = max(abs(real(spec_obs))); end
if data_scale < eps, data_scale = 1; end

spec_norm = spec_obs / data_scale;
sigma_n   = sigma_n  / data_scale;
y_obs_re  = real(spec_norm(mask));
if fit_imag
    y_obs = [y_obs_re; imag(spec_norm(mask))];
else
    y_obs = y_obs_re;
end

%% ---- Amplitude bounds ---------------------------------------------------
if isempty(amp_bounds)
    amp_lb = -0.1 * ones(M,1);
    amp_ub =  Inf(M,1);
else
    amp_lb = amp_bounds(1,:)';
    amp_ub = amp_bounds(2,:)';
end
for k = 1:numel(nonneg_mets)
    ii = find(strcmp(basis.metabolites, nonneg_mets{k}));
    if isempty(ii), error('nonneg_mets: unknown metabolite ''%s''.', nonneg_mets{k}); end
    amp_lb(ii) = max(amp_lb(ii), 0);
end
% Baseline coefficients: unconstrained (one block per fitted channel)
n_bl_tot = n_bl * (1 + fit_imag);
bl_lb = -Inf(n_bl_tot, 1);
bl_ub =  Inf(n_bl_tot, 1);

% Combined bounds for lsqlin: [amps; bl_coeffs]
lin_lb = [amp_lb; bl_lb];
lin_ub = [amp_ub; bl_ub];

%% ---- Ratio constraints (metabolite amplitudes only) ---------------------
A_ineq = zeros(0, M+n_bl_tot);
b_ineq = zeros(0, 1);

if use_constraints || ~isempty(ratio_cons_in)
    mets  = basis.metabolites;
    midx  = @(name) find(strcmp(mets,name));
    if ~isempty(ratio_cons_in)
        if size(ratio_cons_in,2) ~= 4
            error('ratio_constraints must be an N x 4 cell {num, den, lo, hi}.');
        end
        constraints = ratio_cons_in;
    else
        constraints = {
            'NADplus',  'NADH',     5.0, 8.0;
            'PCr',      'bATP',     1.5, 3.0;
            'Pi_intra', 'Pi_extra', 1.0, 5.0;
        };
    end
    for k = 1:size(constraints,1)
        i = midx(constraints{k,1}); j = midx(constraints{k,2});
        lo = constraints{k,3};      hi = constraints{k,4};
        if isempty(i)||isempty(j)
            error('ratio_constraints: unknown metabolite in row %d.', k);
        end
        if ~isempty(ratio_cons_in)
            fprintf('  Ratio constraint: %g <= %s/%s <= %g\n', lo, constraints{k,1}, constraints{k,2}, hi);
        end
        if isfinite(hi)
            row_hi = zeros(1, M+n_bl_tot); row_hi(i)= 1; row_hi(j)=-hi;
            A_ineq = [A_ineq; row_hi]; b_ineq = [b_ineq; 0]; %#ok<AGROW>
        end
        if isfinite(lo)
            row_lo = zeros(1, M+n_bl_tot); row_lo(i)=-1; row_lo(j)= lo;
            A_ineq = [A_ineq; row_lo]; b_ineq = [b_ineq; 0]; %#ok<AGROW>
        end
    end
end

%% ---- Per-group frequency shift setup -----------------------------------
% Groups that get individual chemical shift offsets (bounded to +-fsh_grp_bound).
% PCr+ATP use the global freq_shift (well-determined, used as anchor).
% Other groups can deviate slightly to account for pH, temperature, B0 effects.
%
% Each group offset is RELATIVE to the global freq_shift.
% Final shift for group g: fsh_global + fsh_g
%
% (bound now set by the 'fsh_grp_bound' option; default 0.3 ppm)

% Map each metabolite to a group index (0 = global only, no extra param)
% Groups: 1=PME(PE/PC), 2=PDE(GPC/GPE), 3=Pi, 4=NAD
grp_names  = {'PME','PDE','Pi','NAD','aATP','gATP'};
grp_mets   = {{'PE','PC'}, {'GPC','GPE'}, {'Pi_intra','Pi_extra'}, {'NADH','NADplus'}, {'aATP'}, {'gATP'}};
n_grp      = numel(grp_names);

% Build metabolite->group lookup (0 = global shift only)
met_grp = zeros(M,1);
for g = 1:n_grp
    mets_in_grp = grp_mets{g};          % cell array of metabolite name strings
    for k = 1:numel(mets_in_grp)
        idx = find(strcmp(basis.metabolites, mets_in_grp{k}));
        if ~isempty(idx), met_grp(idx) = g; end
    end
end

% Groups 5-6 (aATP, gATP) are OFF unless asked for, so that old calls
% (scalar or 4-element bound, [4 x 2] prior) behave exactly as before.
if isscalar(fsh_grp_bound_in)
    fsh_grp_bound = [fsh_grp_bound_in * ones(4,1); 0; 0];
elseif numel(fsh_grp_bound_in) == 4
    fsh_grp_bound = [fsh_grp_bound_in(:); 0; 0];
elseif numel(fsh_grp_bound_in) == n_grp
    fsh_grp_bound = fsh_grp_bound_in(:);
else
    error('fsh_grp_bound must be scalar or have 4 or %d elements (%s).', ...
          n_grp, strjoin(grp_names,', '));
end

% Soft priors on group offsets (Hz)
grp_prior = NaN(n_grp, 2);
if ~isempty(grp_prior_in)
    if isequal(size(grp_prior_in), [4 2])
        grp_prior(1:4,:) = grp_prior_in;
    elseif isequal(size(grp_prior_in), [n_grp 2])
        grp_prior = grp_prior_in;
    else
        error('fsh_grp_prior_hz must be [4 x 2] or [%d x 2] = [mean sd] in Hz (NaN row = none).', n_grp);
    end
end
if ~isempty(atp_prior_in)
    if ~isequal(size(atp_prior_in), [2 2])
        error('atp_shift_prior_hz must be [2 x 2] = [mean sd] in Hz for {aATP; gATP}.');
    end
    for q = 1:2
        if all(~isnan(atp_prior_in(q,:))) && atp_prior_in(q,2) > 0
            grp_prior(4+q,:)   = atp_prior_in(q,:);
            fsh_grp_bound(4+q) = (abs(atp_prior_in(q,1)) + 4*atp_prior_in(q,2)) / B0;
        end
    end
end
if nad_rel_aatp
    if fsh_grp_bound(5) <= 0
        error('nad_rel_aatp needs an aATP shift (atp_shift_prior_hz row 1, or fsh_grp_bound(5) > 0).');
    end
    fprintf('  NAD offset measured from the aATP shift (nad_rel_aatp)\n');
end
grp_free = fsh_grp_bound > 0;          % groups whose offset is fitted
n_free   = sum(grp_free);

fprintf('  Per-group freq shifts: %s  (bounds %s ppm)\n', strjoin(grp_names,' | '), ...
        mat2str(fsh_grp_bound'));
has_gprior = ~isnan(grp_prior(:,1)) & ~isnan(grp_prior(:,2)) & grp_prior(:,2) > 0;
for g = find(has_gprior(:))'
    if ~grp_free(g)
        warning('Prior on group %s ignored: its offset is fixed (fsh_grp_bound = 0).', grp_names{g});
        has_gprior(g) = false;
    elseif abs(grp_prior(g,1))/B0 > fsh_grp_bound(g)
        warning('Prior mean for group %s (%.1f Hz) is outside its bound (+-%.1f Hz).', ...
                grp_names{g}, grp_prior(g,1), fsh_grp_bound(g)*B0);
    else
        fprintf('  Prior on %s offset: %.1f +- %.1f Hz\n', grp_names{g}, grp_prior(g,1), grp_prior(g,2));
    end
end

%% ---- Auto-estimate lw_inhom from PCr -----------------------------------
pcr_T2_lw    = 6.0;
pcr_reg_mask = abs(ppm) < 4;
ppm_pcr      = ppm(pcr_reg_mask);
spec_pcr     = real(spec_norm(pcr_reg_mask));
f_pcr_hz     = ppm_pcr * B0;
[A0, pk_idx] = max(spec_pcr);
f0_0         = f_pcr_hz(pk_idx);
lorentz_model = @(p,f) p(1).*(p(3)/2)^2./((f-p(2)).^2+(p(3)/2)^2);
try
    p_fit = lsqcurvefit(lorentz_model,[A0,f0_0,20],f_pcr_hz,spec_pcr, ...
        [0,f0_0-300,1],[Inf,f0_0+300,300], ...
        optimoptions('lsqcurvefit','Display','off'));
    lw_inhom_est = max(2, abs(p_fit(3)) - pcr_T2_lw);
    if lw_inhom_est > 150, lw_inhom_est = 20; end
    fprintf('  PCr Lorentzian fit: FWHM=%.1f Hz -> lw_inhom_est=%.1f Hz\n', ...
            abs(p_fit(3)), lw_inhom_est);
catch
    lw_inhom_est = 20;
    fprintf('  PCr fit failed, fallback lw_inhom=%.1f Hz\n', lw_inhom_est);
end

%% ---- Auto zero-order phase estimate ------------------------------------
% Angle of the PCr peak itself (+-0.15 ppm around the largest point within
% +-0.5 ppm), after removing the local offset measured 0.5-1.5 ppm away, so a
% broad baseline under PCr cannot flip the estimate. The model multiplies the
% basis by exp(+1i*phi0), hence phi0 = angle(peak).
[~, ipk]  = max(abs(spec_norm) .* (abs(ppm) < 0.5));
pk_mask   = abs(ppm - ppm(ipk)) < 0.15;
ring      = abs(ppm - ppm(ipk)) >= 0.5 & abs(ppm - ppm(ipk)) < 1.5;
loc_off   = median(real(spec_norm(ring))) + 1i*median(imag(spec_norm(ring)));
phi0_est  = angle(sum(spec_norm(pk_mask) - loc_off));

%% ---- Initial nonlinear params ------------------------------------------
if isempty(x0_in)
    x0 = [phi0_est; 0; lw_inhom_est; zeros(n_grp,1)];
else
    grp_fsh0 = zeros(n_grp,1);
    if isfield(x0_in,'fsh_groups')
        g0 = x0_in.fsh_groups(:); grp_fsh0(1:min(numel(g0),n_grp)) = g0(1:min(numel(g0),n_grp));
        if nad_rel_aatp, grp_fsh0(4) = grp_fsh0(4) - grp_fsh0(5); end   % effective -> relative
    end
    x0 = [x0_in.phi0; x0_in.freq_shift; x0_in.lw_inhom_hz; grp_fsh0];
end

% Only FREE group offsets enter the nonlinear vector
x0 = [x0(1:3); min(max(x0(3+find(grp_free)), -fsh_grp_bound(grp_free)), fsh_grp_bound(grp_free))];
lb = [-pi; -0.5; lw_inhom_bounds(1); -fsh_grp_bound(grp_free)];
ub = [ pi;  0.5; lw_inhom_bounds(2);  fsh_grp_bound(grp_free)];
i_grp = 3 + (1:n_free);                % positions of free group offsets

% Optional local zero-order phase for a metabolite group (LAST parameter)
lph_idx = [];
for k = 1:numel(lph_mets)
    ii = find(strcmp(basis.metabolites, lph_mets{k}));
    if isempty(ii), error('local_phase_mets: unknown metabolite ''%s''.', lph_mets{k}); end
    lph_idx(end+1) = ii; %#ok<AGROW>
end
has_lph = ~isempty(lph_idx);
i_lph = [];
if has_lph
    lph0 = 0;
    if ~isempty(x0_in) && isfield(x0_in,'local_phase'), lph0 = x0_in.local_phase; end
    x0 = [x0; min(max(lph0, -lph_bound), lph_bound)];
    lb = [lb; -lph_bound];
    ub = [ub;  lph_bound];
    i_lph = numel(x0);
end

% Regional phase / linewidth groups (each one nonlinear parameter + prior)
[pg_names, pg_idx, pg_sd, pg_mu] = parse_groups(phase_groups_in, basis.metabolites, 'phase_groups');
[lg_names, lg_idx, lg_sd, lg_mu] = parse_groups(lw_groups_in,    basis.metabolites, 'lw_groups');
pg_sd = pg_sd * pi/180; pg_mu = pg_mu * pi/180;   % deg -> rad
n_pg = numel(pg_idx); n_lg = numel(lg_idx);
if has_lph && n_pg > 0 && any(ismember(lph_idx, [pg_idx{:}]))
    error('A metabolite is in both local_phase_mets and phase_groups.');
end
i_pg = []; i_lg = [];
for j = 1:n_pg
    x0 = [x0; pg_mu(j)]; lb = [lb; max(pg_mu(j) - 4*pg_sd(j), -pi)]; ub = [ub; min(pg_mu(j) + 4*pg_sd(j), pi)]; %#ok<AGROW>
    i_pg(end+1) = numel(x0); %#ok<AGROW>
    fprintf('  Phase group %-10s : prior %.1f +- %.1f deg (%s)\n', pg_names{j}, pg_mu(j)*180/pi, pg_sd(j)*180/pi, ...
            strjoin(basis.metabolites(pg_idx{j}), ','));
end
for j = 1:n_lg
    x0 = [x0; lg_mu(j)]; lb = [lb; lg_mu(j) - 4*lg_sd(j)]; ub = [ub; lg_mu(j) + 4*lg_sd(j)]; %#ok<AGROW>
    i_lg(end+1) = numel(x0); %#ok<AGROW>
    fprintf('  LW group    %-10s : prior %.1f +- %.1f Hz (%s)\n', lg_names{j}, lg_mu(j), lg_sd(j), ...
            strjoin(basis.metabolites(lg_idx{j}), ','));
end

% Optional t0 (phi1) with a Gaussian prior centred on the given 't0'
fit_t0 = ~isempty(t0_prior_sd);
i_t0 = [];
if fit_t0
    x0 = [x0; t0_fixed];
    lb = [lb; t0_fixed - 4*t0_prior_sd];
    ub = [ub; t0_fixed + 4*t0_prior_sd];
    i_t0 = numel(x0);
    fprintf('  t0 fitted with prior %.1f +- %.1f us\n', t0_fixed*1e6, t0_prior_sd*1e6);
end
i_eta = [];
if fit_eta
    x0 = [x0; eta]; lb = [lb; 0]; ub = [ub; 1];
    i_eta = numel(x0);
    if isempty(eta_prior), fprintf('  eta fitted (0..1), start %.2f\n', eta);
    else, fprintf('  eta fitted with prior %.2f +- %.2f\n', eta_prior); end
end

%% ---- Parameter scaling -------------------------------------------------
scale_fsh = 10;
scale_lw  = 20;
% Per-group shifts use same scale as global fsh
scale_t0  = 1e4;    % t0 in s (~1e-4) -> O(1)
sv   = [1; scale_fsh; 1/scale_lw; scale_fsh*ones(n_free,1); ...
        ones(double(has_lph),1); ones(n_pg,1); (1/scale_lw)*ones(n_lg,1); ...
        scale_t0*ones(double(fit_t0),1); ones(double(fit_eta),1)];
x0_s = x0 .* sv;
lb_s = lb .* sv;
ub_s = ub .* sv;

lsqlin_opts = optimoptions('lsqlin','Display','off','Algorithm','interior-point');

% Baseline curvature penalty (second differences of spline coefficients)
Pen = [];
if ~isempty(bl_curv_sd)
    D2 = diff(eye(n_bl), 2);                       % (n_bl-2) x n_bl
    D2 = D2 / (bl_curv_sd * sigma_n);              % residual in noise units ...
    D2 = D2 * sigma_n;                             % ... expressed like data rows
    Pen = kron(eye(1 + fit_imag), D2);
    fprintf('  Baseline curvature prior: SD %.3g x noise per knot\n', bl_curv_sd);
end

%% ---- VARPRO objective --------------------------------------------------
% Build augmented design matrix [B_metabolites | P_baseline]
% Solve for all linear params (amplitudes + baseline) jointly via lsqlin.
    function f = eff_shift(f)
        % Effective group offsets vs the global shift (NAD rides on aATP if asked)
        if nad_rel_aatp, f(4) = f(4) + f(5); end
    end

    function [res, amps_out, bl_out, B_met] = eval_fit_scaled(xs)
        phi0_x    = xs(1);
        fsh_x     = xs(2) / scale_fsh;
        lw_inh_x  = xs(3) * scale_lw;
        fsh_grp_x = zeros(n_grp,1);          % fixed groups stay at 0
        fsh_grp_x(grp_free) = xs(i_grp) / scale_fsh;
        if fit_t0
            phi1_x = 2*pi * B0 * xs(i_t0) / scale_t0;
        else
            phi1_x = phi1_fixed;
        end

        % Metabolite basis columns (complex) — pass per-group shifts
        lw_off_x = zeros(M,1);
        for jq = 1:n_lg, lw_off_x(lg_idx{jq}) = xs(i_lg(jq)) * scale_lw; end
        if fit_eta, eta_x = xs(i_eta); else, eta_x = eta; end
        B_c   = build_basis_matrix(basis, mask, B0, eta_x, ...
                                   phi0_x, phi1_x, fsh_x, lw_inh_x, ...
                                   met_grp, eff_shift(fsh_grp_x), lw_off_x);
        for jq = 1:n_pg
            B_c(:, pg_idx{jq}) = B_c(:, pg_idx{jq}) * exp(1i*xs(i_pg(jq)));
        end
        if has_lph
            B_c(:, lph_idx) = B_c(:, lph_idx) * exp(1i*xs(i_lph));  % local phase
        end
        B_met = real(B_c);

        % Augmented: metabolites + baseline (+ imaginary channel if requested)
        if fit_imag
            Z = zeros(size(P_w));
            B_aug = [B_met,      P_w, Z; ...
                     imag(B_c),  Z,   P_w];
        else
            B_aug = [B_met, P_w];
        end

        % Regularised baseline: extra rows penalising spline curvature
        if ~isempty(Pen)
            B_aug = [B_aug; zeros(size(Pen,1), M), Pen];
            y_lin = [y_obs; zeros(size(Pen,1),1)];
        else
            y_lin = y_obs;
        end

        % Solve linear params
        c = lsqlin(B_aug, y_lin, A_ineq, b_ineq, [], [], ...
                   lin_lb, lin_ub, [], lsqlin_opts);

        amps_out = c(1:M);
        bl_out   = c(M+1:end);
        res      = B_aug * c - y_lin;          % data rows first, then penalty rows

        % Soft priors, scaled to noise units of the data residual
        pr = zeros(0,1);
        for gq = find(has_gprior(:))'
            pr(end+1,1) = sigma_n * (fsh_grp_x(gq)*B0 - grp_prior(gq,1)) / grp_prior(gq,2); %#ok<AGROW>
        end
        if fit_t0
            pr(end+1,1) = sigma_n * (xs(i_t0)/scale_t0 - t0_fixed) / t0_prior_sd;
        end
        for jq = 1:n_pg
            pr(end+1,1) = sigma_n * (xs(i_pg(jq)) - pg_mu(jq)) / pg_sd(jq); %#ok<AGROW>
        end
        for jq = 1:n_lg
            pr(end+1,1) = sigma_n * (xs(i_lg(jq)) * scale_lw - lg_mu(jq)) / lg_sd(jq); %#ok<AGROW>
        end
        if fit_eta && ~isempty(eta_prior)
            pr(end+1,1) = sigma_n * (xs(i_eta) - eta_prior(1)) / eta_prior(2);
        end
        res = [res; pr];
    end

%% ---- Nonlinear optimisation --------------------------------------------
nl_opts = optimoptions('lsqnonlin', ...
    'Display',                  'off',   ...
    'MaxFunctionEvaluations',   10000,   ...
    'MaxIterations',            1000,    ...
    'FunctionTolerance',        1e-10,   ...
    'StepTolerance',            1e-10,   ...
    'OptimalityTolerance',      1e-10,   ...
    'FiniteDifferenceType',     'central', ...
    'FiniteDifferenceStepSize', 1e-3);

fprintf('  Starting optimisation: phi0=%.3f rad  fsh=%.4f ppm  lw_inhom=%.2f Hz\n', ...
        x0(1), x0(2), x0(3));

[x_opt_s, resnorm_val, ~, exitflag] = ...
    lsqnonlin(@eval_fit_scaled, x0_s, lb_s, ub_s, nl_opts);

x_opt = x_opt_s ./ sv;

%% ---- Recover parameters ------------------------------------------------
phi0_fit    = x_opt(1);
fsh_fit     = x_opt(2);
lw_inh_fit  = x_opt(3);
fsh_grp_fit = zeros(n_grp,1);
fsh_grp_fit(grp_free) = x_opt(i_grp);   % [n_grp x 1], fixed groups = 0
nad_rel_fit = fsh_grp_fit(4);           % NAD parameter (relative to aATP if nad_rel_aatp)
fsh_grp_fit = eff_shift(fsh_grp_fit);   % effective offsets vs the global shift
if fit_t0, t0_fit = x_opt(i_t0); else, t0_fit = t0_fixed; end
phi1_fit = 2*pi * B0 * t0_fit;
if has_lph, lph_fit = x_opt(i_lph); else, lph_fit = 0; end
if fit_eta, eta_fit = x_opt(i_eta); fprintf('  eta fitted: %.3f\n', eta_fit); else, eta_fit = eta; end
pg_fit = x_opt(i_pg); lg_fit = x_opt(i_lg);
met_phase = zeros(M,1);  met_phase(lph_idx) = lph_fit;
for j = 1:n_pg, met_phase(pg_idx{j}) = pg_fit(j); end
met_lwoff = zeros(M,1);
for j = 1:n_lg, met_lwoff(lg_idx{j}) = lg_fit(j); end

[res_final, amps_fit_norm, bl_coeffs_norm, B_fit] = eval_fit_scaled(x_opt_s);
n_obs    = numel(y_obs);
rss_data = sum(res_final(1:n_obs).^2);              % data part only
prior_pen = sum(res_final(n_obs+1:end).^2) / sigma_n^2;   % chi2 units
k_total  = M + n_bl_tot + numel(x0);                % all fitted parameters

amps_fit   = amps_fit_norm * data_scale;
amps_raw   = amps_fit;                       % unclipped solution
bl_coeffs_norm_im = bl_coeffs_norm(n_bl+1:end);   % empty unless fit_imag
bl_coeffs_norm    = bl_coeffs_norm(1:n_bl);
bl_coeffs  = bl_coeffs_norm * data_scale;
if clip_negative
    amps_fit(amps_fit < 0) = 0;
end

lw_total_lines = basis.prior.lw_T2_hz + lw_inh_fit;

fprintf('  Amplitudes:\n');
for m = 1:M
    fprintf('    %-12s : %.4e\n', basis.metabolites{m}, amps_fit(m));
end
fprintf('  phi0=%.3f rad  fsh=%.4f ppm  lw_inhom=%.2f Hz  exitflag=%d\n', ...
        phi0_fit, fsh_fit, lw_inh_fit, exitflag);
if has_lph
    fprintf('  Local phase (%s): %+.4f rad (%+.1f deg)\n', ...
            strjoin(lph_mets, ','), lph_fit, lph_fit*180/pi);
end
for j = 1:n_pg
    fprintf('  Phase group %-10s : %+.1f deg\n', pg_names{j}, pg_fit(j)*180/pi);
end
for j = 1:n_lg
    fprintf('  LW group    %-10s : %+.2f Hz\n', lg_names{j}, lg_fit(j));
end
if fit_t0
    fprintf('  t0 fitted: %.1f us (prior %.1f +- %.1f us) -> phi1=%.4f rad/ppm\n', ...
            t0_fit*1e6, t0_fixed*1e6, t0_prior_sd*1e6, phi1_fit);
end
if prior_pen > 0
    fprintf('  Prior penalty (chi2 units): %.2f\n', prior_pen);
end
fprintf('  Per-group freq offsets (ppm):\n');
for g = 1:n_grp
    if g > 4 && fsh_grp_bound(g) == 0, continue; end
    fprintf('    %-6s : %+.4f ppm  (%+.1f Hz)\n', ...
            grp_names{g}, fsh_grp_fit(g), fsh_grp_fit(g)*B0);
end
if nad_rel_aatp
    fprintf('    NAD - aATP : %+.1f Hz (fitted relative offset)\n', nad_rel_fit*B0);
end

%% ---- Reconstruct -------------------------------------------------------
phase_vec = exp(1i*(phi0_fit + phi1_fit.*ppm));
spec_c    = zeros(N,1);
spec_each = zeros(N, M);           % real part of each metabolite's fitted spectrum
for m = 1:M
    g   = met_grp(m);
    fsh_m = fsh_fit + (g>0)*fsh_grp_fit(max(g,1));   % global + group offset
    idx = strcmp(basis.prior.metabolite, basis.metabolites{m});
    pk  = basis.prior(idx,:);
    fid_m = zeros(N,1);
    for r = 1:height(pk)
        cs_Hz    = (pk.cs_ppm(r) + fsh_m) * B0;
        lw_total = max(pk.lw_T2_hz(r) + lw_inh_fit + met_lwoff(m), 0.5);
        fid_m    = fid_m + amps_fit(m) * make_component(t, cs_Hz, lw_total, pk.rel_amp(r), eta_fit);
    end
    spec_m = fftshift(fft(fid_m)) .* phase_vec * exp(1i*met_phase(m));
    spec_each(:, m) = real(spec_m);
    spec_c = spec_c + spec_m;
end
spec_metabs_full = real(spec_c);
baseline_full    = real(P_full * bl_coeffs);
if fit_imag
    baseline_imag_full = real(P_full * bl_coeffs_norm_im) * data_scale;
else
    baseline_imag_full = [];
end
spec_fit_full    = spec_metabs_full + baseline_full;
residual_full    = real(spec_obs) - spec_fit_full;

%% ---- Metrics (metabolite columns only — baseline not in FIM) -----------
metrics = compute_metrics(y_obs_re - real(P_w * bl_coeffs_norm), ...
    B_fit, amps_fit_norm, sigma_n, ...
    spec_metabs_full/data_scale, ppm, ppm_range, ...
    basis.metabolites, rss_data, n_fit, M, k_total, n_obs);

metrics.CRLB_abs      = metrics.CRLB_abs     * data_scale;
metrics.posterior_std = metrics.posterior_std * data_scale;
metrics.sigma_noise   = sigma_n * data_scale;
metrics.CRLB_pct      = 100 * metrics.CRLB_abs ./ (abs(amps_fit) + eps);
metrics.table.Amplitude = amps_fit(:);
metrics.table.CRLB_abs  = metrics.CRLB_abs(:);
metrics.table.CRLB_pct  = metrics.CRLB_pct(:);

%% ---- Pack result -------------------------------------------------------
result.amplitudes      = amps_fit;
result.amplitudes_raw  = amps_raw;          % unclipped (may be negative)
result.fsh_grp_bound   = fsh_grp_bound;
result.fit_imag        = fit_imag;
result.local_phase     = lph_fit;           % rad (0 if option off)
result.local_phase_mets = lph_mets;
result.fit_t0          = fit_t0;
result.phase_group_names = pg_names;
result.phase_group_deg = pg_fit(:)' * 180/pi;
result.lw_group_names  = lg_names;
result.lw_group_hz     = lg_fit(:)';
result.met_phase       = met_phase;         % rad, all metabolites
result.met_lw_offset   = met_lwoff;
result.eta             = eta_fit;
result.baseline_curv_sd = bl_curv_sd;         % Hz,  all metabolites
result.rss_data        = rss_data * data_scale^2;   % data residual only
result.prior_penalty   = prior_pen;
result.n_params        = k_total;
result.baseline_imag_full = baseline_imag_full;
result.phi0            = phi0_fit;
result.t0              = t0_fit;
result.phi1            = phi1_fit;
result.freq_shift      = fsh_fit;
result.fsh_groups      = fsh_grp_fit;   % [n_grp x 1] per-group offsets (ppm)
result.fsh_group_names = grp_names;
result.nad_rel_aatp    = nad_rel_aatp;
result.nad_rel_aatp_hz = nad_rel_fit * B0;   % = NAD offset vs global if nad_rel_aatp is false
result.lw_inhom_hz     = lw_inh_fit;
result.lw_total_hz     = lw_total_lines;
result.baseline_coeffs = bl_coeffs;
result.baseline_full   = baseline_full;
result.spec_fit        = spec_fit_full;
result.spec_metabs     = spec_metabs_full;
result.spec_each       = spec_each;          % [N x M], one column per metabolite
result.spec_obs        = real(spec_obs);
result.residual        = residual_full;
result.ppm             = ppm;
result.metabolites     = basis.metabolites;
result.metrics         = metrics;
result.exitflag        = exitflag;
result.resnorm         = resnorm_val;
result.sigma_noise     = sigma_n * data_scale;
result.data_scale      = data_scale;
result.baseline_order  = K;   % stored as n_knots for display

if doPlot, plot_result(result, basis); end
end

%% =========================================================================
function B = build_basis_matrix(basis, mask, B0, eta, phi0, phi1, fsh, lw_inh, met_grp, fsh_grp, lw_off)
% Returns [n_fit x M] matrix of COMPLEX metabolite basis spectra
% (callers take real() / imag() as needed).
% met_grp : [M x 1] group index per metabolite (0 = global shift only)
% fsh_grp : [n_grp x 1] per-group frequency offsets (ppm), relative to fsh
if nargin < 9, met_grp = zeros(numel(basis.metabolites),1); end
if nargin < 10, fsh_grp = zeros(max(met_grp(:)),1); end
if nargin < 11, lw_off = zeros(numel(basis.metabolites),1); end

N  = numel(basis.t);
t  = basis.t(:);
M  = numel(basis.metabolites);
pv = exp(1i*(phi0 + phi1.*basis.ppm(:)));
B  = zeros(sum(mask), M);
for m = 1:M
    g     = met_grp(m);
    fsh_m = fsh + (g>0 && ~isempty(fsh_grp)) * fsh_grp(max(g,1));
    idx   = strcmp(basis.prior.metabolite, basis.metabolites{m});
    pk    = basis.prior(idx,:);
    fid_m = zeros(N,1);
    for r = 1:height(pk)
        cs_Hz    = (pk.cs_ppm(r) + fsh_m) * B0;
        lw_total = max(pk.lw_T2_hz(r) + lw_inh + lw_off(m), 0.5);
        fid_m    = fid_m + make_component(t, cs_Hz, lw_total, pk.rel_amp(r), eta);
    end
    sp_m   = fftshift(fft(fid_m)) .* pv;
    B(:,m) = sp_m(mask);
end
end

%% =========================================================================
function [names, idx, sd, mu] = parse_groups(G, mets, optname)
% G: N x 3 cell {name, {metabolites}, prior_sd} or N x 4 with prior mean.
names = {}; idx = {}; sd = zeros(0,1); mu = zeros(0,1);
if isempty(G), return; end
if ~any(size(G,2) == [3 4])
    error('%s must be an N x 3 (or N x 4) cell: {name, {metabolites}, prior_sd [, prior_mean]}.', optname);
end
used = [];
for j = 1:size(G,1)
    ml = G{j,2}; if ischar(ml), ml = {ml}; end
    ii = zeros(1, numel(ml));
    for k = 1:numel(ml)
        f = find(strcmp(mets, ml{k}));
        if isempty(f), error('%s: unknown metabolite ''%s''.', optname, ml{k}); end
        ii(k) = f;
    end
    if any(ismember(ii, used)), error('%s: a metabolite appears in two groups.', optname); end
    if ~(isscalar(G{j,3}) && G{j,3} > 0), error('%s: prior SD must be a positive scalar.', optname); end
    m0 = 0;
    if size(G,2) == 4 && ~isempty(G{j,4}), m0 = G{j,4}; end
    used = [used, ii]; %#ok<AGROW>
    names{end+1} = G{j,1}; idx{end+1} = ii; sd(end+1,1) = G{j,3}; mu(end+1,1) = m0; %#ok<AGROW>
end
end

%% =========================================================================
function s = make_component(t, cs_Hz, lw_hz, amp, eta)
R2s   = pi * lw_hz;
sig_t = sqrt(2*log(2)) / (pi*lw_hz + eps);
lor   = exp(-R2s .* t);
gau   = exp(-t.^2 ./ (2*sig_t^2));
s     = amp .* exp(1i*2*pi*cs_Hz.*t) .* ((1-eta)*lor + eta*gau);
end

%% =========================================================================
function metrics = compute_metrics(y_bl_corrected, B, amps, sigma_n, ...
    spec_fit, ppm, ppm_range, metabolites, resnorm, n_fit, M, k_total, n_obs)
% y_bl_corrected: observed minus fitted baseline — for FIM on metabolites only

% CRLB = sqrt(diag(inv(FIM))) with FIM = B'B/sigma^2, i.e. sigma*sqrt(diag(inv(B'B))).
% (Earlier versions multiplied by sigma a second time -> CRLBs ~1/sigma_norm
%  too small, typically 50-100x.) Metabolite columns only: baseline and
%  nonlinear parameters are not included, so these are still optimistic.
FIM_inv  = pinv(B'*B);
crlb_abs = sigma_n * sqrt(max(diag(FIM_inv), 0));
crlb_pct = 100 * crlb_abs ./ (abs(amps) + eps);

mask  = ppm >= ppm_range(1) & ppm <= ppm_range(2);
res_v = y_bl_corrected - real(spec_fit(mask));
residual_norm = norm(res_v) / (norm(y_bl_corrected) + eps);

fit_SNR = zeros(M,1);
fit_FNR = zeros(M,1);
for m = 1:M
    fit_SNR(m) = max(abs(B(:,m))) * abs(amps(m)) / (sigma_n + eps);
    fit_FNR(m) = sum(abs(B(:,m))) * abs(amps(m)) / (sigma_n*sqrt(n_fit) + eps);
end

% All fitted parameters; n_obs = number of fitted data points (x2 if imag)
% Standard Gaussian-likelihood forms (additive constants dropped)
n_params = k_total;
AIC = n_obs*log(resnorm/n_obs) + 2*n_params;
BIC = n_obs*log(resnorm/n_obs) + n_params*log(n_obs);

metrics.CRLB_abs      = crlb_abs;
metrics.CRLB_pct      = crlb_pct;
metrics.residual_norm = residual_norm;
metrics.fit_SNR       = fit_SNR;
metrics.fit_FNR       = fit_FNR;
metrics.AIC           = AIC;
metrics.BIC           = BIC;
metrics.posterior_std = crlb_abs;
metrics.sigma_noise   = sigma_n;
metrics.metabolites   = metabolites;
metrics.table = table(metabolites(:), amps(:), crlb_abs(:), crlb_pct(:), ...
    fit_SNR(:), fit_FNR(:), ...
    'VariableNames',{'Metabolite','Amplitude','CRLB_abs','CRLB_pct','fitSNR','fitFNR'});
end

%% =========================================================================
function plot_result(result, basis)
ppm    = result.ppm;
N      = numel(ppm);
M      = numel(result.metabolites);
colors = lines(M);
t      = basis.t(:);
B0     = basis.B0_MHz;
eta    = basis.eta; if isfield(result, 'eta'), eta = result.eta; end
pv     = exp(1i*(result.phi0 + result.phi1.*ppm));

LB_hz = 10;
k_ap  = (0:N-1).';
apod  = exp(-pi * LB_hz .* k_ap / basis.fs);
apfun = @(s) real(fftshift(fft(ifft(ifftshift(s)) .* apod)));

% All quantities are already in frequency domain — mild Gaussian smooth for display
gsmooth = @(s) smoothdata(s, 'gaussian', 3);

obs_d  = gsmooth(result.spec_obs);
met_d  = gsmooth(result.spec_metabs);
bl_d   = gsmooth(result.baseline_full);
fit_d  = gsmooth(result.spec_fit);
res_d  = gsmooth(result.residual);

figure('Name','31P Fit','Color','w','Position',[60 60 1150 720]);
ax1 = subplot(3,1,1:2); hold on;
plot(ppm, obs_d, 'Color',[0.65 0.65 0.65],'LineWidth',0.8,'DisplayName','Observed');
plot(ppm, fit_d, 'r','LineWidth',1.8,'DisplayName','Total fit');
plot(ppm, bl_d,  'Color',[0.85 0.55 0.1],'LineStyle','-.','LineWidth',1.4, ...
     'DisplayName',sprintf('Baseline (order %d)', result.baseline_order));

for m = 1:M
    g     = 0;
    if isfield(result,'fsh_groups') && isfield(result,'fsh_group_names')
        % rebuild met_grp locally for plotting
        grp_mets_plot = {{'PE','PC'}, {'GPC','GPE'}, {'Pi_intra','Pi_extra'}, {'NADH','NADplus'}, {'aATP'}, {'gATP'}};
        for gg = 1:numel(grp_mets_plot)
            if any(strcmp(grp_mets_plot{gg}, result.metabolites{m})), g = gg; break; end
        end
    end
    fsh_m = result.freq_shift + (g>0)*result.fsh_groups(max(g,1));
    idx   = strcmp(basis.prior.metabolite, result.metabolites{m});
    pk    = basis.prior(idx,:);
    fid_m = zeros(numel(t),1);
    for r = 1:height(pk)
        cs_Hz    = (pk.cs_ppm(r) + fsh_m) * B0;
        lw_total = pk.lw_T2_hz(r) + result.lw_inhom_hz;
        if isfield(result,'met_lw_offset'), lw_total = max(lw_total + result.met_lw_offset(m), 0.5); end
        fid_m    = fid_m + make_component(t, cs_Hz, lw_total, pk.rel_amp(r), eta);
    end
    pv_m = pv;
    if isfield(result,'met_phase')
        pv_m = pv * exp(1i*result.met_phase(m));
    elseif isfield(result,'local_phase_mets') && any(strcmp(result.local_phase_mets, result.metabolites{m}))
        pv_m = pv * exp(1i*result.local_phase);
    end
    sp_m = real(fftshift(fft(result.amplitudes(m)*fid_m)) .* pv_m);
    plot(ppm, gsmooth(sp_m), '--', 'Color',colors(m,:), 'LineWidth',0.9, ...
         'DisplayName', sprintf('%s (SNR=%.1f)', result.metabolites{m}, ...
                                result.metrics.fit_SNR(m)));
end

set(ax1,'XDir','reverse');
xlabel('ppm'); ylabel('Intensity');
title(sprintf('31P Fit  |  exitflag=%d  |  LW_{inhom}=%.2f Hz  |  t_0=%.3f ms  |  BL knots=%d', ...
      result.exitflag, result.lw_inhom_hz, result.t0*1e3, result.baseline_order));
legend('Location','northwest','FontSize',7,'NumColumns',3);
xlim([min(ppm) max(ppm)]); box off; grid on;

ax2 = subplot(3,1,3); hold on;
plot(ppm, res_d, 'Color',[0.3 0.3 0.3],'LineWidth',0.8);
yline(0,'k--','LineWidth',0.6);
set(ax2,'XDir','reverse');
xlabel('ppm'); ylabel('Residual');
title(sprintf('Residual (norm=%.4f)', result.metrics.residual_norm));
xlim([min(ppm) max(ppm)]); box off; grid on;
linkaxes([ax1 ax2],'x');

% Metrics table
figure('Name','Fit Metrics','Color','w','Position',[1200 60 800 440]);
T   = result.metrics.table;
fmt = @(v) arrayfun(@(x) sprintf('%.4e',x), v, 'UniformOutput',false);
dat = [T.Metabolite, fmt(T.Amplitude), fmt(T.CRLB_abs), ...
       fmt(T.CRLB_pct), fmt(T.fitSNR), fmt(T.fitFNR)];
uitable('Data',dat, ...
    'ColumnName',{'Metabolite','Amplitude','CRLB_abs','CRLB_%','fitSNR','fitFNR'}, ...
    'Units','normalized','Position',[0.02 0.18 0.96 0.76],'FontSize',10);
annotation('textbox',[0.02 0.01 0.96 0.12], ...
    'String',sprintf('AIC=%.1f  BIC=%.1f  Residual norm=%.4f  LW_{inhom}=%.2f Hz  BL knots=%d', ...
        result.metrics.AIC, result.metrics.BIC, result.metrics.residual_norm, ...
        result.lw_inhom_hz, result.baseline_order), ...
    'EdgeColor','none','HorizontalAlignment','center','FontSize',10);
end

%% =========================================================================
function B = bspline_basis(x, t, n)
% BSPLINE_BASIS  Evaluate n cubic (order-4) B-spline basis functions at x.
%
% Uses the de Boor recursion (Cox-de Boor formula):
%   N_{i,1}(x) = 1 if t(i) <= x < t(i+1), else 0
%   N_{i,k}(x) = (x-t(i))/(t(i+k-1)-t(i)) * N_{i,k-1}(x)
%              + (t(i+k)-x)/(t(i+k)-t(i+1)) * N_{i+1,k-1}(x)
%
% x : [nx x 1] evaluation points
% t : augmented knot vector (length n+4 for cubic)
% n : number of basis functions
% B : [nx x n] basis matrix  (each column is locally supported)

x  = x(:);
nx = numel(x);
p  = 4;   % cubic = order 4

% Build order-1 basis (indicator functions)
N = zeros(nx, numel(t)-1);
for i = 1:numel(t)-1
    if i < numel(t)-1
        N(:,i) = (x >= t(i)) & (x < t(i+1));
    else
        % Last interval: close the right endpoint
        N(:,i) = (x >= t(i)) & (x <= t(i+1));
    end
end

% Recurse up to order p
for k = 2:p
    N_new = zeros(nx, numel(t)-k);
    for i = 1:numel(t)-k
        d1 = t(i+k-1) - t(i);
        d2 = t(i+k)   - t(i+1);
        c1 = 0; c2 = 0;
        if d1 > 0, c1 = (x - t(i))   / d1 .* N(:,i);   end
        if d2 > 0, c2 = (t(i+k) - x) / d2 .* N(:,i+1); end
        N_new(:,i) = c1 + c2;
    end
    N = N_new;
end

B = N(:, 1:n);

% Normalise columns to [0,1] range for numerical stability in lsqlin
col_max = max(abs(B), [], 1);
col_max(col_max < eps) = 1;
B = B ./ col_max;
end