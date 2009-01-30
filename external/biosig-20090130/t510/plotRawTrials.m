function plotRawTrials(r)

trials = r{1}.trials;
triallen = r{1}.triallen;
st = r{1}.st;
cue = r{1}.cue;
class = r{1}.class;
plim = r{1}.plim;
name = r{1}.name;
fs = r{1}.fs;

% Plot raw trials
fig = figure('Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);
set(gcf, 'DefaultAxesFontSize', 8);

for k = 1:trials
    subplot(trials, 1, k);
    plot(st(1, (k-1)*triallen+1:k*triallen));
    set(gca, 'XTick', 0:fs:triallen, 'TickLength', [0.001 0.001]);
    set(gca, 'XTickLabel', []);
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');
    ylabel(num2str(k));
    v = axis;

    if ~isempty(cue)
        line([cue*fs cue*fs], [v(3) v(4)], 'Color', 'k');
    end;
    if ~isempty(plim)
        axis([v(1) v(2) plim]);
    end;

end;

set(gca, 'XTickLabel', 0:triallen/fs);

[nrows, ncols] = size(class);
if (ncols == 1)
    class = class';
end;

% Edit plot
annotation('textbox', [0 0 0.3 0.1], 'LineStyle', 'none', 'Interpreter', ...
           'none', 'String', {[name ', ' date], ['fs = ' num2str(fs) ' Hz'], ...
           ['trials = ' num2str(trials) ', length = ' num2str(triallen/fs) ...
           ' s'], ['class = ' num2str(class)]});

hndl = findobj('parent', fig, 'type', 'axes');
for a = 1:length(hndl)
    set(findobj('parent', hndl(a)), 'ButtonDownFcn', 'zoomRawTrials');
end;