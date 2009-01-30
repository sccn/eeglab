function plotHR(r)

cl = r{1}.cl;
cu = r{1}.cu;
hrt = r{1}.hrt;
triallen = r{1}.triallen;
trials = r{1}.trials;
alpha = r{1}.alpha;
cue = r{1}.cue;
name = r{1}.name;
class = r{1}.class;
fs = r{1}.fs;

fig = figure('Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);
set(gcf, 'DefaultAxesFontSize', 8);

plot(mean(hrt, 1), 'LineWidth', 2);
hold on;
if (~isempty(alpha))
    plot(cl);
    plot(cu);
end;
line([0 triallen], [mean(mean(hrt)) mean(mean(hrt))], 'Color', 'k');

set(gca, 'XTick', 0:fs:triallen);
set(gca, 'XTickLabel', 0:triallen/fs);

if (~isempty(cue))
    v = axis;
    line([cue*fs cue*fs], [v(3) v(4)], 'Color', 'k');
end;

[nrows, ncols] = size(class);
if (ncols == 1)
    class = class';
end;

% Edit plot
annotation('textbox', [0 0 0.3 0.1], 'LineStyle', 'none', 'Interpreter', ...
           'none', 'String', {[name ', ' date], ['fs = ' num2str(fs) ' Hz'], ...
           ['trials = ' num2str(trials) ', length = ' num2str(triallen/fs) ...
           ' s'], ['class = ' num2str(class)]});