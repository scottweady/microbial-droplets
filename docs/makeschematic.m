% MAKESCHEMATIC  Light and dark versions of the Fig. 1 schematic for the docs page.
%
% Rasterizes manuscript/fig1.pdf with Ghostscript, crops the empty margins and
% writes docs/assets/fig1-light.png. The schematic mixes vector art with an
% embedded CMYK raster, so the dark version is made from the rendered image:
% in CIELAB the lightness is reversed onto the page's dark-theme range
% (white -> background #1f2424, black -> text #ececec) while the hue (a*, b*)
% is kept, so lines become light, faces dark, and antialiasing is preserved.
% Writes docs/assets/fig1-dark.png.
%
% Run from the repository root (needs Ghostscript, `gs`):
%   matlab -batch "run('docs/makeschematic.m')"

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
src = fullfile(root, 'manuscript', 'fig1.pdf');
out = fullfile(here, 'assets');

tmp = [tempname '.png'];
gs = '/usr/local/bin/gs';
if ~isfile(gs), gs = 'gs'; end
status = system(sprintf(['%s -q -dNOPAUSE -dBATCH -sDEVICE=png16m -r220 ' ...
  '-dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sOutputFile="%s" "%s"'], gs, tmp, src));
assert(status == 0, 'Ghostscript failed')
I = im2double(imread(tmp));
delete(tmp)

% Crop to the drawn content plus a small margin
ink = any(I < 0.97, 3);
rows = find(any(ink, 2)); cols = find(any(ink, 1));
pad = 16;
I = I(max(rows(1) - pad, 1) : min(rows(end) + pad, end), ...
      max(cols(1) - pad, 1) : min(cols(end) + pad, end), :);
% Transparent surround: the white region connected to the image border,
% widened by a few pixels to take in the antialiased edges of the outlines.
% Inside the box (fluid, droplet) the image stays opaque.
white = all(I > 0.985, 3);
outside = imfill(~white, 'holes');               % box and everything drawn
outside = ~outside;                              % white surround only
band = imdilate(outside, strel('disk', 4));

% Dark theme: reverse lightness, keep hue, lift chroma slightly
bgRGB = [31 36 36] / 255;                        % --bg
bg = rgb2lab(bgRGB);
fg = rgb2lab([236 236 236] / 255);               % --fg
L = rgb2lab(I);
s = 1 - L(:, :, 1) / 100;                        % 0 on white, 1 on black
D = zeros(size(L));
D(:, :, 1) = bg(1) + s * (fg(1) - bg(1));
D(:, :, 2) = 1.3 * L(:, :, 2) + (1 - s) * bg(2);
D(:, :, 3) = 1.3 * L(:, :, 3) + (1 - s) * bg(3);
D = min(max(lab2rgb(D), 0), 1);

[Il, Al] = colorToAlpha(I, [1 1 1], band);
[Id, Ad] = colorToAlpha(D, bgRGB, band);
imwrite(Il, fullfile(out, 'fig1-light.png'), 'Alpha', Al);
imwrite(Id, fullfile(out, 'fig1-dark.png'), 'Alpha', Ad);

fprintf('wrote docs/assets/fig1-light.png and docs/assets/fig1-dark.png\n');


function [C, A] = colorToAlpha(I, b, mask)
  % Remove the background colour b from the pixels in mask (as GIMP's
  % colour-to-alpha): the smallest alpha such that compositing the result
  % over b reproduces I exactly. Pixels outside mask stay opaque.
  b = reshape(b, 1, 1, 3);
  up = (I - b) ./ max(1 - b, eps);                % channels above the background
  dn = (b - I) ./ max(b, eps);                    % channels below it
  a = max(max(up, dn), [], 3);
  a = min(max(a, 0), 1);
  A = ones(size(a));
  A(mask) = a(mask);
  C = I;
  aa = max(A, 1e-6);
  C = b + (I - b) ./ aa;
  C(repmat(A < 1e-6, 1, 1, 3)) = 0;
  C = min(max(C, 0), 1);
end
