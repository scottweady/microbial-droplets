% MAKEINTERACTIVE  Data for the interactive figures on the docs page.
%
% Evaluates the axisymmetric solution on the unit disk (beta = 1, kappa = 0)
% and writes
%   docs/assets/surface-data.js  window.SURFACE = { r, ug, ub }
%       radial surface velocity outside the droplet, 1 <= r <= 3, for growth
%       (Ra = 0) and for buoyancy at unit beta*Ra (Fig. 3(b), Fig. 8(b)).
%   docs/assets/bulk-data.js     window.BULK = { x, z, ug, wg, ub, wb, c, ice, matter }
%       bulk velocity in the plane y = 0 for growth and for unit beta*Ra
%       buoyancy (Fig. 4), and the nutrient concentration c (Fig. 2).
%   docs/assets/kappa-data.js    window.KAPPA = { H, kappa }   (Fig. 9)
%
% The bulk kernels are those of figs/scripts/3DFlow_plot.jl. The total
% solution at any beta*Ra is the linear combination  u = u_g + beta*Ra * u_b.
% (.js rather than .json so that the page also works when opened from disk.)
%
% Run from the repository root:  matlab -batch "run('docs/makeinteractive.m')"

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
assets = fullfile(here, 'assets');
if ~isfolder(assets), mkdir(assets); end

lam10 = 4 / pi;                               % lambda_1^0

% ── Quadrature on the unit disk ─────────────────────────────────────────────
% r = sin(phi) clusters nodes at the edge, where p ~ sqrt(1 - r^2).
nphi = 240; nth = 720;
[phi, wphi] = gaussLegendre(nphi, 0, pi/2);
th = 2 * pi * (0 : nth - 1)' / nth;
[PHI, TH] = ndgrid(phi, th);
W = (sin(PHI) .* cos(PHI) .* wphi) * (2 * pi / nth);   % r dr dtheta
R = sin(PHI);
yx = R .* cos(TH); yy = R .* sin(TH);
pg = lam10 * cos(PHI);                                 % growth pressure
pb = -(lam10 / 96) * cos(PHI) .* (1 + (4/3) * cos(PHI).^2);  % buoyancy pressure
sig = 2 * ones(size(PHI));                             % density, beta = 1
yx = yx(:); yy = yy(:); W = W(:); pg = pg(:); pb = pb(:); sig = sig(:);

% ── Surface velocity outside the droplet (z = 0) ────────────────────────────
rext = [1, linspace(1.01, 3, 240)];
ugs = zeros(size(rext)); ubs = zeros(size(rext));
for n = 2 : numel(rext)
  X1 = rext(n) - yx; X2 = -yy;
  d2 = X1.^2 + X2.^2; d = sqrt(d2);
  Gg = X1 ./ (4 * pi * d.^3);
  Gb = X1 / (128 * pi) .* (1 - 2 * log(d));
  ugs(n) = sum(W .* pg .* Gg);
  ubs(n) = sum(W .* (pb .* Gg + sig .* Gb));
end
ugs(1) = 0.5; ubs(1) = 0;                              % continuity at r = 1
fprintf('surface: u_g(1.01) = %.4f, max 96 u_b = %.4f at r = %.3f\n', ...
        ugs(2), 96 * max(ubs), rext(find(ubs == max(ubs), 1)));

S = struct('r', rnd(rext), 'ug', rnd(ugs), 'ub', rnd(ubs));
writeJS(fullfile(assets, 'surface-data.js'), 'SURFACE', S, ...
  'radial surface velocity outside the unit droplet, growth and unit beta*Ra buoyancy');

% ── Bulk velocity and concentration in the plane y = 0 ──────────────────────
nxh = 65; nz = 65;
xh = linspace(0, 2, nxh);                              % x >= 0, mirrored below
z = linspace(-2, -0.02, nz);
UG = zeros(nz, nxh); WG = UG; UB = UG; WB = UG; C = UG;
for i = 1 : nz
  zi = z(i);
  for j = 1 : nxh
    X1 = xh(j) - yx; X2 = -yy;
    s2 = X1.^2 + X2.^2;
    rho = sqrt(s2 + zi^2);
    rmz = rho - zi;                                    % > 0 since z < 0
    gp = (1 ./ rho.^3 - 3 * zi^2 ./ rho.^5) / (4 * pi);
    Ggx = X1 .* gp;
    Ggz = zi * (s2 - 2 * zi^2) ./ (4 * pi * rho.^5);
    Gbx = X1 / (128 * pi) .* (1 + 4 * zi ./ rho - 2 * log(rmz)) ...
          - 6 * zi * X1 ./ (128 * pi * rmz);
    Gbz = zi / (32 * pi) * (zi ./ rho + log(rmz));
    UG(i, j) = sum(W .* pg .* Ggx);
    WG(i, j) = sum(W .* pg .* Ggz);
    UB(i, j) = sum(W .* (pb .* Ggx + sig .* Gbx));
    WB(i, j) = sum(W .* (pb .* Ggz + sig .* Gbz));
    C(i, j)  = 1 - sum(W .* sig ./ (4 * pi * rho));
  end
end
% Mirror to x < 0: u is odd in x, w and c are even
x = [-fliplr(xh(2:end)), xh];
mir = @(A, s) [s * fliplr(A(:, 2:end)), A];
UG = mir(UG, -1); UB = mir(UB, -1); WG = mir(WG, 1); WB = mir(WB, 1); C = mir(C, 1);
fprintf('bulk: u_g(0.5, z=%.2f) = %.4f (surface value 0.25)\n', z(end), UG(end, x == 0.5));

flat = @(A) rnd(reshape(A.', 1, []));                  % row-major, x fastest
B = struct('x', rnd(x), 'z', rnd(z), ...
           'ug', flat(UG), 'wg', flat(WG), 'ub', flat(UB), 'wb', flat(WB), ...
           'c', flat(C), ...
           'ice', round(cmocean('ice', 64), 4), ...
           'matter', round(flipud(cmocean('matter', 64)), 4));
writeJS(fullfile(assets, 'bulk-data.js'), 'BULK', B, ...
  'bulk velocity (growth, unit beta*Ra buoyancy) and concentration, plane y = 0');

% ── Matching constant kappa against domain depth |H| (Fig. 9) ───────────────
K = load(fullfile(root, 'data', 'kappa_H.dat'));
keep = K(:, 1) >= 0.5 & K(:, 1) <= 8 & isfinite(K(:, 2));
K = K(keep, :);
K = K(1 : 2 : end, :);
writeJS(fullfile(assets, 'kappa-data.js'), 'KAPPA', ...
  struct('H', rnd(K(:, 1)'), 'kappa', rnd(K(:, 2)')), ...
  'matching constant kappa against effective depth |H|');


% ── Helpers ─────────────────────────────────────────────────────────────────
function v = rnd(v)
  v = round(v, 5, 'significant');
end

function writeJS(file, name, s, what)
  fid = fopen(file, 'w');
  fprintf(fid, '// Generated by docs/makeinteractive.m: %s.\n', what);
  fprintf(fid, 'window.%s = %s;\n', name, jsonencode(s));
  fclose(fid);
end

function [x, w] = gaussLegendre(n, a, b)
  % Golub-Welsch
  k = 1 : n - 1;
  beta = k ./ sqrt(4 * k.^2 - 1);
  [V, D] = eig(diag(beta, 1) + diag(beta, -1));
  [x, i] = sort(diag(D));
  w = 2 * V(1, i)'.^2;
  x = (b - a) / 2 * x + (a + b) / 2;
  w = (b - a) / 2 * w;
end
