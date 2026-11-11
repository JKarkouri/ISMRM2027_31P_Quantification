function [b, info] = make_basis_ren2020(basis, varargin)
% MAKE_BASIS_REN2020  Rebuild the NAD(H) and UDP-sugar entries of a 31P basis
% from Ren, Malloy & Sherry, MRM 2020;84:2338 (Table 1 shifts, J couplings),
% with every component normalised to its number of 31P nuclei.
%
%   b = make_basis_ren2020(basis)
%   b = make_basis_ren2020(basis, 'udp_mode','separate', 'shift_set','IR')
%
% Normalisation: one 31P nucleus = the summed rel_amp of PCr (normally 1).
% NADH (singlet, 2 equivalent P) and NAD+ (AB quartet) then both sum to 2,
% each UDP-sugar (two weakly coupled doublets, Palpha + Pbeta) to 2, so an
% amplitude ratio is a molar ratio, and amplitude/PCr amplitude = conc/[PCr].
%
% Options
%   'shift_set'   'conventional' (default; fit with aATP in the model, as
%                 ours) or 'IR' (their inversion-edited fit). Table 1.
%   'udp_mode'    'combined' (default): one metabolite 'UDPG' = the four
%                 UDP-sugars with fixed proportions from Table 2.
%                 'separate': four metabolites UDPGal, UDPGlc, UDPGalNAc,
%                 UDPGlcNAc (replace 'UDPG').  'keep': leave UDPG as is
%                 (only renormalised).
%   'nadh_shift'  override NADH shift (ppm)                     default []
%   'normalise_all' true -> also rescale every other metabolite so its
%                 rel_amp sums to its number of 31P (PCr, Pi, PE, PC, GPC,
%                 GPE and each ATP group = 1). Amplitude / PCr amplitude is
%                 then a concentration ratio for every metabolite. Does not
%                 change the fit, only the amplitude scale.   default false
%
% The NADH linewidth is NOT changed here; Ren et al. found NADH ~1.7x broader
% than NAD+; use 'lw_groups' in the fitter, e.g. {'NADH',{'NADH'},8}.
%
% info: table of each metabolite's summed rel_amp before and after.

p = inputParser;
addParameter(p, 'shift_set', 'conventional', @(s) any(strcmpi(s, {'conventional','IR'})));
addParameter(p, 'udp_mode',  'combined',     @(s) any(strcmpi(s, {'combined','separate','keep'})));
addParameter(p, 'nadh_shift', [], @(x) isempty(x) || isscalar(x));
addParameter(p, 'normalise_all', false, @islogical);
parse(p, varargin{:});
o = p.Results;

B0 = basis.B0_MHz;
T  = basis.prior;
mets_before = basis.metabolites(:)';
sum_before = cellfun(@(m) sum(T.rel_amp(strcmp(T.metabolite, m))), mets_before);

perP = sum(T.rel_amp(strcmp(T.metabolite, 'PCr')));
if isempty(perP) || perP <= 0, error('PCr not found in basis.prior.'); end

% ---- Ren et al. 2020, Table 1 (ppm, PCr = 0) and J (Hz) -------------------
if strcmpi(o.shift_set, 'IR')
    sh.NADH = -8.092; sh.NADp = [-8.132 -8.467];
    sh.Gal = [-8.047 -9.554]; sh.Glc = [-8.165 -9.697];
    sh.GalNAc = [-8.215 -9.785]; sh.GlcNAc = [-8.265 -9.961];
    frac = [0.048 0.089 0.079 0.116];          % Table 2, IR multi-doublet (mM)
else
    sh.NADH = -8.111; sh.NADp = [-8.131 -8.467];
    sh.Gal = [-8.022 -9.543]; sh.Glc = [-8.162 -9.699];
    sh.GalNAc = [-8.218 -9.782]; sh.GlcNAc = [-8.265 -9.963];
    frac = [0.047 0.091 0.080 0.117];          % Table 2, conventional multi-doublet
end
if ~isempty(o.nadh_shift), sh.NADH = o.nadh_shift; end
J.NADp = 20.14; J.UDP = [20.3 20.3 20.5 20.5];
udp_names = {'UDPGal','UDPGlc','UDPGalNAc','UDPGlcNAc'};
udp_sh = {sh.Gal, sh.Glc, sh.GalNAc, sh.GlcNAc};

% ---- template rows (keep any extra columns, e.g. lw_T2_hz) ---------------
tpl = @(name) T(find(strcmp(T.metabolite, name), 1), :);
if ~any(strcmp(T.metabolite, 'NADH')) || ~any(strcmp(T.metabolite, 'NADplus'))
    error('basis must contain NADH and NADplus.');
end
rowNADH = tpl('NADH'); rowNADp = tpl('NADplus');
hasUDPG = any(strcmp(T.metabolite, 'UDPG'));
if hasUDPG, rowUDP = tpl('UDPG'); else, rowUDP = rowNADp; end

% ---- new rows ------------------------------------------------------------
new = T([], :);
new = [new; mkrows(rowNADH, 'NADH', sh.NADH, 2*perP)];
[pos, amp] = ab_lines(sh.NADp(1), sh.NADp(2), J.NADp, B0);
new = [new; mkrows(rowNADp, 'NADplus', pos, amp * perP)];      % sums to 2*perP

mets = basis.metabolites(:)';
switch lower(o.udp_mode)
    case 'combined'
        w = frac / sum(frac);
        for k = 1:4
            [pos, amp] = ab_lines(udp_sh{k}(1), udp_sh{k}(2), J.UDP(k), B0);
            new = [new; mkrows(rowUDP, 'UDPG', pos, amp * perP * w(k))]; %#ok<AGROW>
        end
        if ~hasUDPG, mets{end+1} = 'UDPG'; end
        drop = {'NADH','NADplus','UDPG'};
    case 'separate'
        for k = 1:4
            [pos, amp] = ab_lines(udp_sh{k}(1), udp_sh{k}(2), J.UDP(k), B0);
            new = [new; mkrows(rowUDP, udp_names{k}, pos, amp * perP)]; %#ok<AGROW>
        end
        iU = find(strcmp(mets, 'UDPG'));
        if isempty(iU)
            mets = [mets, udp_names];
        else
            mets = [mets(1:iU-1), udp_names, mets(iU+1:end)];
        end
        drop = {'NADH','NADplus','UDPG'};
    otherwise   % keep: renormalise UDPG only
        if hasUDPG
            r = T(strcmp(T.metabolite, 'UDPG'), :);
            r.rel_amp = r.rel_amp / sum(r.rel_amp) * 2 * perP;
            new = [new; r];
        end
        drop = {'NADH','NADplus','UDPG'};
end

keep = ~ismember(T.metabolite, drop);
b = basis;
b.prior = [T(keep, :); new];
b.metabolites = mets;

if o.normalise_all
    nP1 = {'PCr','Pi_intra','Pi_extra','PE','PC','GPC','GPE','gATP','aATP','bATP'};
    for k = 1:numel(nP1)
        ii = strcmp(b.prior.metabolite, nP1{k});
        if any(ii), b.prior.rel_amp(ii) = b.prior.rel_amp(ii) / sum(b.prior.rel_amp(ii)) * perP; end
    end
end

% ---- report --------------------------------------------------------------
sum_after = cellfun(@(m) sum(b.prior.rel_amp(strcmp(b.prior.metabolite, m))), mets);
all_m = unique([mets_before, mets], 'stable');
sb = nan(numel(all_m),1); sa = sb;
for k = 1:numel(all_m)
    i1 = strcmp(mets_before, all_m{k}); if any(i1), sb(k) = sum_before(i1); end
    i2 = strcmp(mets, all_m{k});        if any(i2), sa(k) = sum_after(i2);  end
end
info = table(all_m(:), sb, sa, 'VariableNames', {'metabolite','sum_rel_amp_before','sum_rel_amp_after'});
fprintf('\nmake_basis_ren2020 (%s shifts, UDP %s), 1 P = %.3g:\n', o.shift_set, o.udp_mode, perP);
disp(info);
nadp = b.prior(strcmp(b.prior.metabolite,'NADplus'), :);
fprintf('NAD+ lines (ppm): %s, rel_amp %s\n', mat2str(nadp.cs_ppm', 4), mat2str(nadp.rel_amp', 3));
fprintf('NADH line (ppm): %.3f\n', sh.NADH);
end

%% ---------------------------------------------------------------------------
function [pos_ppm, amp] = ab_lines(dA, dB, J, B0)
% Four lines of an AB system (two coupled 31P), total intensity 2.
nuA = dA*B0; nuB = dB*B0; c = (nuA+nuB)/2;
D = sqrt((nuA-nuB)^2 + J^2);
pos_hz = c + [ (D+J)/2, (D-J)/2, -(D-J)/2, -(D+J)/2 ];
amp    = [1-J/D, 1+J/D, 1+J/D, 1-J/D] / 2;
pos_ppm = pos_hz / B0;
end

function R = mkrows(tplrow, name, pos, amp)
% Copy a template row once per line, set name / shift / rel_amp.
n = numel(pos);
R = repmat(tplrow, n, 1);
if iscell(R.metabolite)
    R.metabolite = repmat({name}, n, 1);
else
    R.metabolite = repmat(string(name), n, 1);
end
R.cs_ppm  = pos(:);
R.rel_amp = amp(:);
end
