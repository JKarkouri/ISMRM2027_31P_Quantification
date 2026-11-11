function A = amares_fit_spectrum(spec, b, pk, info, varargin)
% AMARES_FIT_SPECTRUM  Fit one spectrum with AMARES (OXSA) and return amounts
% per metabolite plus the curves for the figures.
%
%   [pk, info] = make_amares_pk(bNew);
%   A = amares_fit_spectrum(spec_p, bNew, pk, info);
%
% spec : fftshifted spectrum on the basis ppm axis (PCr at 0 ppm, carrier at 0
%        ppm, nominal delay already removed), as used for the other fitters.
% Output fields:
%   A.mets        metabolites (basis order, excluded ones absent)
%   A.area        summed line amplitudes per metabolite (FID amplitude units)
%   A.per_mol     area / 31P per molecule (NAD+, NADH, UDP(G): 2)
%   A.crlb_pct    CRLB of each metabolite's amplitude, % (first line)
%   A.fit, A.comp fitted spectrum (metabolites + A.bg) and one spectrum
%                 per metabolite (data frame, i.e. with AMARES' fitted phases),
%                 A.data = spec
%   A.bg          baseline: the broad peaks (make_amares_pk 'broad_ppm') plus,
%                 with skip_points, the part of the data in the skipped points
%   A.res, A.status  raw AMARES outputs
% Options: 'beginTime' (0 s), 'MaxIter' (500), 'quiet' (true),
%   'skip_points' (0): leave out the first K FID points (AMARES then models the
%   FID from t = beginTime + K*dt). This is the usual AMARES way of removing
%   broad components (baseline): they have decayed after ~1 ms. The part left
%   out is returned as A.bg, the spectrum of the residual over those K points,
%   so data = fit + bg + residual as for the fitters with a baseline.

p = inputParser;
addParameter(p, 'beginTime', 0);
addParameter(p, 'MaxIter', 500);
addParameter(p, 'quiet', true);
addParameter(p, 'skip_points', 0);
parse(p, varargin{:}); o = p.Results;

N = numel(spec); t = b.t(:) - b.t(1); dt = t(2) - t(1);
K = o.skip_points;
ex.samples = N - K; ex.dwellTime = dt; ex.imagingFrequency = b.B0_MHz;
ex.timeAxis = t(1:N-K); ex.ppmAxis = b.ppm(:); ex.beginTime = o.beginTime + K*dt;
ex.offset = 0; ex.expectedOffset = 0;
fid = ifft(ifftshift(double(spec(:))));
[res, st, ~, crb] = AMARES.amaresFit(fid(K+1:end), ex, pk, false, 'quiet', o.quiet, 'MaxIter', o.MaxIter);

lm = info.lineMet(:)';
if isfield(info, 'broad'), isB = ismember(lm, info.broad); else, isB = false(size(lm)); end
mets = unique(lm(~isB), 'stable');
if numel(res.amplitude) ~= numel(lm)
    error('AMARES returned %d lines, the prior knowledge describes %d.', numel(res.amplitude), numel(lm));
end
[~, fids] = AMARES.makeModelFid(res, t + o.beginTime, b.B0_MHz);
A.mets = mets; A.area = zeros(1, numel(mets)); A.per_mol = A.area; A.crlb_pct = A.area;
A.comp = zeros(N, numel(mets));
for m = 1:numel(mets)
    j = strcmp(lm, mets{m});
    A.area(m) = sum(res.amplitude(j));
    A.per_mol(m) = A.area(m) / info.nP(mets{m});
    jj = find(j, 1);
    A.crlb_pct(m) = 100 * crb.amplitude(jj) / max(res.amplitude(jj), eps);
    A.comp(:, m) = fftshift(fft(sum(fids(:, j), 2)));
end
A.data = spec(:);
broad = fftshift(fft(sum(fids(:, isB), 2)));     % broad (baseline) peaks, if any
rf = fid - sum(fids, 2);                         % residual FID over all points
A.bg = broad + fftshift(fft([rf(1:K); zeros(N-K, 1)]));  % baseline: broad peaks + what the skipped points carried
A.fit = sum(A.comp, 2) + A.bg;                   % drawn over the data: metabolites + baseline
                                                 % (data - fit = residual of the fitted points only)
A.res = res; A.status = st;
end
