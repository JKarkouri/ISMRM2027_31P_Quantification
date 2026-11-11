function [pk, info] = make_amares_pk(b, varargin)
% MAKE_AMARES_PK  AMARES (OXSA) prior knowledge built from a basis struct, so
% AMARES fits the same metabolites, chemical shifts, J-couplings and relative
% line amplitudes as the in-house fitter, LCModel and FSL-MRS.
%
%   [pk, info] = make_amares_pk(bNew);                       % corrected (Ren) basis
%   [pk, info] = make_amares_pk(bNew, 'write', 'PK_7T_Brain_Ren2020_31P.m');
%   [pk, info] = make_amares_pk(basis_real);                 % original basis
%
% How the basis lines become AMARES peaks
%   - Each metabolite's lines are split into "sites" (gaps > site_gap_ppm), e.g.
%     the two phosphates of UDP(G). Each site is one AMARES (multiplet) peak.
%   - 1 line  -> singlet; 2 lines -> doublet; 3 lines -> triplet (outer/centre
%     ratio from the basis); 4 lines (NAD+ AB quartet) -> quartet with the
%     outer/inner ratio from the basis and the one spacing that best matches the
%     unequal AB spacings; > 4 lines (UDP(G): several sugars) -> doublet at the
%     weighted centres. info.max_err_hz gives the worst line-position error.
%   - AMARES multiplet convention: chemShift = FIRST (lowest-ppm) line,
%     lines at chemShift + k*chemShiftDelta.
%   - Sites of the same metabolite share amplitude (if their 31P count is
%     equal) and linewidth (G_amplitude, G_linewidth).
%   - Phases ('phase_mode'): 'nad_aatp' (default) = each peak has its own
%     zero-order phase, except that the sites of one molecule share one, and
%     NAD+, NADH and the UDP-sugars share the alpha-ATP phase (small peaks on
%     the alpha-ATP shoulder, which a free phase would let absorb baseline);
%     'all' = one phase for all peaks; 'free' = every peak its own phase.
%     (The old 'link_phases' true/false = 'all'/'free'.)
%   - info.lineMet{i} is the metabolite of the i-th fitted line (AMARES expands
%     multiplets in row order); info.nP(met) the 31P per molecule, so
%     per-molecule amount = sum of the metabolite's line amplitudes / nP.
%
% Options: 'shift_tol_ppm' (0.1), 'small_tol_ppm' (0.05, NAD/UDP), 'lw_bounds'
% ([1 200] Hz), 'small_lw_bounds' ([10 40] Hz), 'lw_extra' (10 Hz added to the
% basis width as start), 'phase_mode' ('nad_aatp'), 'exclude' ({'MP_broad'}),
% 'broad_ppm' ([] = none; e.g. [0 -2.5] adds two broad baseline peaks),
% 'broad_lw' ([150 3000] Hz), 'broad_tol_ppm' (4),
% 'site_gap_ppm' (0.6), 'write' ('' or a .m file name for the PK function).

p = inputParser;
addParameter(p, 'shift_tol_ppm', 0.1);
addParameter(p, 'small_tol_ppm', 0.05);
addParameter(p, 'small_mets', {'NADH','NADplus','UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'});
addParameter(p, 'lw_bounds', [1 200]);
addParameter(p, 'small_lw_bounds', [10 40]);
addParameter(p, 'lw_extra', 10);
addParameter(p, 'phase_mode', 'nad_aatp');      % 'nad_aatp' | 'all' | 'free'
addParameter(p, 'link_phases', []);             % old option: true = 'all', false = 'free'
addParameter(p, 'phase_to_aatp', {'NADH','NADplus','UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'});
addParameter(p, 'exclude', {'MP_broad'});
addParameter(p, 'broad_ppm', []);            % broad (baseline) peaks: centres in ppm, e.g. [0 -2.5]
addParameter(p, 'broad_lw', [150 3000]);     % Hz, bounds of their linewidths
addParameter(p, 'broad_tol_ppm', 4);         % ppm, bounds of their shifts
addParameter(p, 'lineshape', 'lorentz');     % 'voigt' = OXSA Lorentz x Gauss, one shared Gaussian width
addParameter(p, 'sigma_bounds', [0 30]);     % Hz, Gaussian sigma (voigt only)
addParameter(p, 'sigma0', 3);                % Hz, its start value
addParameter(p, 'site_gap_ppm', 0.6);
addParameter(p, 'write', '');
addParameter(p, 'nad_quartet', 'two_doublets');    % or 'equal_spacing' (classic AMARES quartet, approximate AB positions)
addParameter(p, 'lw_link', {{'UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'}, {'NADH','NADplus'}});   % metabolites sharing one linewidth
parse(p, varargin{:}); o = p.Results;

mets = b.metabolites(:)';
mets = mets(~ismember(mets, o.exclude));
twoP = {'NADH','NADplus','UDPG','UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
info.nP = containers.Map();
R = struct('names', {}, 'tag', {}, 'met', {}, 'first', {}, 'mult', {}, 'd', {}, 'r', {}, ...
           'lw0', {}, 'lwb', {}, 'tol', {}, 'w', {}, 'err', {}, 'split', {});
for m = 1:numel(mets)
    msk = strcmp(b.prior.metabolite, mets{m});
    cs = b.prior.cs_ppm(msk); ra = b.prior.rel_amp(msk); lt = b.prior.lw_T2_hz(msk);
    [pos, o2] = sort(cs(:)); w = ra(o2); w = w(:); lw = lt(o2);
    if any(strcmp(mets{m}, twoP)), info.nP(mets{m}) = 2; else, info.nP(mets{m}) = 1; end
    cut = [0; find(diff(pos) > o.site_gap_ppm); numel(pos)];
    small = any(strcmp(mets{m}, o.small_mets));
    for s = 1:numel(cut) - 1
        ii = cut(s)+1:cut(s+1);
        if numel(ii) == 4 && strcmp(o.nad_quartet, 'two_doublets')
            parts = {ii([2 3]), ii([1 4])};         % inner and outer doublet, independent amplitudes
        else
            parts = {ii};
        end
        for q = 1:numel(parts)
            jj = parts{q};
            [first, mult, d, r, model_pos] = site_model(pos(jj), w(jj));
            k = numel(R) + 1;
            tag = mets{m};
            if numel(cut) > 2, tag = sprintf('%s_s%d', tag, s); end
            if numel(parts) > 1, tag = sprintf('%s_%s', tag, char('i' + (q-1)*('o'-'i'))); end
            if isempty(mult), names = tag; else, names = arrayfun(@(j) sprintf('%s_%d', tag, j), 1:numel(mult), 'UniformOutput', false); end
            R(k).names = names; R(k).tag = tag; R(k).met = mets{m}; R(k).first = first; R(k).mult = mult; R(k).d = d; R(k).r = r;
            R(k).lw0 = mean(lw(jj)) + o.lw_extra;
            if small, R(k).lwb = o.small_lw_bounds; R(k).tol = o.small_tol_ppm;
            else, R(k).lwb = o.lw_bounds; R(k).tol = o.shift_tol_ppm; end
            R(k).lw0 = min(max(R(k).lw0, R(k).lwb(1) + 0.5), R(k).lwb(2) - 0.5);
            R(k).w = sum(w(jj)); R(k).split = numel(parts) > 1;
            R(k).err = max(abs(nearest_err(pos(jj), model_pos))) * b.B0_MHz;
        end
    end
end

% groups: rows of the same metabolite share linewidth (and amplitude if equal 31P)
iPCr = find(strcmp({R.met}, 'PCr'), 1);
if isempty(iPCr), iPCr = 1; end
rowsBound = cell(numel(R), 7); rowsIV = cell(numel(R), 5); rowsPK = cell(numel(R), 9);
info.lineMet = {};
for k = 1:numel(R)
    same = find(strcmp({R.met}, R(k).met));
    lead = same(1);
    gLw = []; gAmp = [];
    lwgrp = R(k).met;                              % linewidth group: the metabolite, or its lw_link set
    for g = 1:numel(o.lw_link)
        if any(strcmp(o.lw_link{g}, R(k).met)), lwgrp = o.lw_link{g}; end
    end
    lwrows = find(ismember({R.met}, lwgrp));
    if numel(lwrows) > 1, gLw = lwrows(1); end
    if numel(same) > 1 && ~R(k).split            % sites of one molecule: same 31P count -> same amplitude
        if abs(R(k).w - R(lead).w) < 0.1*R(lead).w, gAmp = lead;
        else, warning('%s: sites have different 31P counts; amplitudes left independent.', R(k).met); end
    end
    gPh = [];                                     % set below when phases are linked
    rowsBound(k,:) = {R(k).names, R(k).first + [-1 1]*R(k).tol, R(k).lwb, [0 Inf], [0 360], [], []};
    rowsIV(k,:)    = {R(k).names, R(k).first, R(k).lw0, 1, 0};
    rowsPK(k,:)    = {R(k).names, R(k).mult, R(k).d, R(k).r, gLw, gAmp, gPh, [], double(k == iPCr)};
    nl = max(1, numel(R(k).mult));
    info.lineMet(end+1:end+nl) = {R(k).met};
end
if ~isempty(o.link_phases)
    if o.link_phases, o.phase_mode = 'all'; else, o.phase_mode = 'free'; end
end
switch o.phase_mode
    case 'all'                                     % all phases in ONE group led by row 1
        rowsPK(:, 7) = {1};
    case 'nad_aatp'
        % groups: the sites of each molecule; NAD/UDP rows joined to the alpha-ATP rows.
        % AMARES: the group value is the index of the group's first (lowest) row.
        iA = find(strcmp({R.met}, 'aATP'));
        key = {R.met};
        if ~isempty(iA), key(ismember(key, o.phase_to_aatp)) = {'aATP'}; end
        for k = 1:numel(R)
            g = find(strcmp(key, key{k}));
            if numel(g) > 1, rowsPK{k, 7} = min(g); end
        end
    case 'free'
    otherwise
        error('phase_mode must be ''nad_aatp'', ''all'' or ''free''.');
end
info.phase_mode = o.phase_mode;

% broad peaks: AMARES' way of fitting a broad baseline (hump, first-point
% distortion). Singlets with wide linewidths and their own free phase (a
% phase of 180 deg gives a negative dip); not linked to anything; reported
% as baseline (A.bg in amares_fit_spectrum), not as metabolites.
info.broad = {};
for j = 1:numel(o.broad_ppm)
    nmB = sprintf('BROAD%d', j);
    rowsBound(end+1, :) = {nmB, o.broad_ppm(j) + [-1 1]*o.broad_tol_ppm, o.broad_lw, [0 Inf], [0 360], [], []}; %#ok<AGROW>
    rowsIV(end+1, :)    = {nmB, o.broad_ppm(j), sqrt(prod(o.broad_lw)), 1, 0}; %#ok<AGROW>
    rowsPK(end+1, :)    = {nmB, [], [], [], [], [], [], [], 0}; %#ok<AGROW>
    info.lineMet{end+1} = nmB; info.nP(nmB) = 1; info.broad{end+1} = nmB;
end
fields.Bounds = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase', 'chemShiftDelta', 'amplitudeRatio'};
fields.IV = {'peakName', 'chemShift', 'linewidth', 'amplitude', 'phase'};
fields.PK = {'peakName', 'multiplet', 'chemShiftDelta', 'amplitudeRatio', 'G_linewidth', 'G_amplitude', ...
             'G_phase', 'G_chemShiftDelta', 'refPeak'};
if strcmpi(o.lineshape, 'voigt')
    % OXSA Voigt: each line x exp(-2 pi^2 sigma^2 t^2). One Gaussian width shared by
    % all metabolite peaks (G_sigma -> row 1), like the single eta of the in-house
    % fitter; broad baseline peaks stay Lorentzian (sigma fixed at 0).
    nR = numel(R); nAll = size(rowsBound, 1);
    rowsBound(:, 8) = [repmat({o.sigma_bounds}, nR, 1); repmat({[]}, nAll - nR, 1)];
    rowsIV(:, 6)    = [repmat({o.sigma0}, nR, 1); repmat({0}, nAll - nR, 1)];
    rowsPK(:, 10)   = [repmat({1}, nR, 1); repmat({[]}, nAll - nR, 1)];
    fields.Bounds{end+1} = 'sigma'; fields.IV{end+1} = 'sigma'; fields.PK{end+1} = 'G_sigma';
end
info.lineshape = lower(o.lineshape);
values.boundsCellArray = rowsBound; values.IVCellArray = rowsIV; values.PKCellArray = rowsPK;
info.rows = R; info.max_err_hz = max([R.err]);

fprintf('AMARES prior knowledge from the basis (%d peaks, %d lines):\n', numel(R), numel(info.lineMet));
fprintf('  %-12s %-10s %8s %8s %8s %8s\n', 'peak', 'type', 'first', 'delta', 'ratio', 'err(Hz)');
for k = 1:numel(R)
    typ = {'singlet','doublet','triplet','quartet'}; t = typ{min(max(1, numel(R(k).mult)), 4)};
    fprintf('  %-12s %-10s %8.3f %8.3f %8.3f %8.1f\n', R(k).tag, t, ...
        R(k).first, empty0(R(k).d), empty0(R(k).r), R(k).err);
end

if ~isempty(o.write), write_pk_file(o.write, fields, values); end
if ~isempty(which('AMARES.priorKnowledge.preparePriorKnowledge'))
    pk = AMARES.priorKnowledge.preparePriorKnowledge(fields, values);
else
    warning('OXSA (AMARES) is not on the path: returning the raw fields/values.');
    pk = struct('fields', fields, 'values', values);
end
end

%% -----------------------------------------------------------------------
function [first, mult, d, r, mp] = site_model(pos, w)
n = numel(pos);
c = sum(pos .* w) / sum(w);
switch n
    case 1
        first = pos; mult = []; d = []; r = []; mp = pos;
    case 2
        first = pos(1); mult = [0 0]; d = pos(2) - pos(1); r = []; mp = pos;
    case 3
        d = mean(diff(pos)); first = c - d; mult = [1 0 1]; r = mean(w([1 3])) / w(2);
        mp = first + (0:2)'*d;
    case 4
        a = (pos(3) - pos(2))/2; bb = (pos(4) - pos(1))/2;
        d = (a + 3*bb)/5;                       % least-squares single spacing for an AB quartet
        first = c - 1.5*d; mult = [1 0 0 1]; r = mean(w([1 4])) / mean(w([2 3]));
        mp = first + (0:3)'*d;
    otherwise                                   % several overlapping doublets (UDP-sugars)
        lo = pos <= c; hi = ~lo;
        cl = sum(pos(lo).*w(lo))/sum(w(lo)); cu = sum(pos(hi).*w(hi))/sum(w(hi));
        first = cl; mult = [0 0]; d = cu - cl; r = []; mp = [cl; cu];
end
end

function e = nearest_err(pos, mp)
e = arrayfun(@(x) min(abs(mp - x)), pos);
end

function v = empty0(x)
if isempty(x), v = 0; else, v = x; end
end

function write_pk_file(fname, fields, values)
[~, fn] = fileparts(fname);
fid = fopen(fname, 'w');
fprintf(fid, 'function outStruct = %s()\n', fn);
fprintf(fid, '%% AMARES prior knowledge generated by make_amares_pk from the basis struct (%s).\n', datestr(now, 'yyyy-mm-dd')); %#ok<TNOW1,DATST>
fprintf(fid, '%% Multiplets: chemShift = first (lowest-ppm) line; lines at chemShift + k*chemShiftDelta.\n\n');
blocks = {'Bounds', 'boundsCellArray'; 'IV', 'IVCellArray'; 'PK', 'PKCellArray'};
for k = 1:3
    fprintf(fid, 'fields.%s = {%s};\n', blocks{k,1}, strjoin(cellfun(@(s) ['''' s ''''], fields.(blocks{k,1}), 'UniformOutput', false), ', '));
    fprintf(fid, 'values.%s = {...\n', blocks{k,2});
    C = values.(blocks{k,2});
    for r = 1:size(C, 1)
        fprintf(fid, '    %s;\n', strjoin(cellfun(@lit, C(r,:), 'UniformOutput', false), ', '));
    end
    fprintf(fid, '    };\n\n');
end
fprintf(fid, 'outStruct = AMARES.priorKnowledge.preparePriorKnowledge(fields, values);\nend\n');
fclose(fid);
fprintf('Wrote %s\n', fname);
end

function s = lit(x)
if ischar(x), s = ['''' x ''''];
elseif iscell(x), s = ['{' strjoin(cellfun(@(c) ['''' c ''''], x, 'UniformOutput', false), ',') '}'];
elseif isempty(x), s = '[]';
elseif isscalar(x), s = sprintf('%.10g', x);
else, s = ['[' strjoin(arrayfun(@(v) sprintf('%.10g', v), x, 'UniformOutput', false), ',') ']'];
end
end
