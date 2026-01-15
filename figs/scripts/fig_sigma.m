
close all
fig = journal_figure([6.5 2.5], 2, 'tex');

% Choose solution: 'g' for growth, 'b' for buoyancy
% sol = 'g';

% Mode numbers
m = 0 : 40;

% Load and set formatting parameters
if strcmp(sol, 'g')
%   sigma = load('../../data/sigma_g.dat');
  sigma = 1/2 - m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4;
  name = '{\it \sigma_g}'; 
  ys = [-4 1]; 
  yt = -4 : 1; 
  loc = 'southwest'; 
  c = [0.4 0.8 0.4];
elseif strcmp(sol, 'b')
  sigma = (1 / 96) * (m .* lambda_PSH(1, 0) .* lambda_PSH(m, m) / 4 - 9 ./ (2 * (m - 1) .* (2 * m + 1)));
  sigma(1:2) = 0;
  name = '{\it \sigma_b}'; 
  ys = [-0.01 0.05]; 
  yt = -0.01 : 0.01 : 0.05; 
  loc = 'northwest'; 
  c = [0.4 0.4 0.8];
end

% Plot
plot(m, sigma, 'o-', 'Color', c, 'DisplayName', name), hold on
plot(m, 0 * m, '--', 'Color', 'k', 'DisplayName', 'none', 'LineWidth', 1, 'HandleVisibility','off')
plot(m, sigma(1) + (sigma(end)-sigma(1)) * (m / m(end)).^0.5, ':', 'Color', [0.2 0.2 0.2], 'DisplayName', '\propto {\it m^{1/2}}')
xlabel('mode, {\it m}');
ylabel(strcat('stability coefficient, ', name))

xlim([0 30])
ylim(ys)
xticks(0 : 5 : 30)
yticks(yt)

legend('location', loc)

% Format
ax = gca;
ax.Units = 'inches';

ax.Position(3) = 0.45 * fig.PaperSize(1);
ax.Position(1) = 0.5 * fig.PaperSize(1) - 0.5 * ax.Position(3);
ax.Position(4) = 0.7 * fig.PaperSize(2);
ax.Position(2) = 0.5 * fig.PaperSize(2) - 0.5 * ax.Position(4);
ax.FontSize = 18;

function val = lambda_PSH(l, m)

  id = abs(m) <= l;
  val = (gamma((l + m + 1) / 2) .* gamma((l - m + 1) / 2)) ./ (gamma((l + m + 2) / 2) .* gamma((l - m + 2) / 2));
  val = val .* id;

end
