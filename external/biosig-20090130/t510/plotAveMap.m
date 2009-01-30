function plotAveMap(r)

plot_index = r{1}.plot_index;
n_rows = r{1}.n_rows;
n_cols = r{1}.n_cols;
fs = r{1}.fs;
name = r{1}.name;
trials = r{1}.trials;
class = r{1}.class;
triallen = r{1}.triallen;
startend = r{1}.startend;
cue = r{1}.cue;

% Plot average maps
fig = figure('Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);
set(gcf, 'DefaultAxesFontSize', 6);

peaksE = [0 0];  % Minimum and maximum of mean (EEG)
peaksM = [0 0];  % Minimum and maximum of mean (EMG)
peaksO = [0 0];  % Minimum and maximum of mean (EOG)


for k = 1:length(plot_index)

    if (r{k}.chantype == 'E')  % All EEG channels have the same scaling
        if (min(r{k}.average) < peaksE(1))
            peaksE(1) = min(r{k}.average);
        end;
        if (max(r{k}.average) > peaksE(2))
            peaksE(2) = max(r{k}.average);
        end;
    end;
    if (r{k}.chantype == 'O')  % All EOG channels have the same scaling
        if (min(r{k}.average) < peaksO(1))
            peaksO(1) = min(r{k}.average);
        end;
        if (max(r{k}.average) > peaksO(2))
            peaksO(2) = max(r{k}.average);
        end;
    end;
    if (r{k}.chantype == 'M')  % All EMG channels have the same scaling
        if (min(r{k}.average) < peaksM(1))
            peaksM(1) = min(r{k}.average);
        end;
        if (max(r{k}.average) > peaksM(2))
            peaksM(2) = max(r{k}.average);
        end;
    end;
end;

for k = 1:length(plot_index)
    subplot(n_rows, n_cols, plot_index(k));

    ax = plotyy(1:triallen, r{k}.average, 1:triallen, r{k}.stdev, 'plot');

    if (r{k}.chantype == 'E')  % All EEG channels have the same scaling
        axis(ax(1), [0 triallen floor(peaksE(1)/10)*10 ceil(peaksE(2)/10)*10]);
        axis(ax(2), [0 triallen 0 100]);
        %set(ax(1), 'YTick', floor(peaksE(1)/10)*10:10:ceil(peaksE(2)/10)*10, 'FontSize', 8);
        %set(ax(1), 'YTickLabel', floor(peaksE(1)/10)*10:10:ceil(peaksE(2)/10)*10, 'FontSize', 8);
        set(ax(2), 'YTick', 0:20:40);
        set(ax(2), 'YTickLabel', 0:20:40);
    elseif (r{k}.chantype == 'M')  % all EMG channels have the same scaling
        axis(ax(1), [0 triallen floor(peaksM(1)/10)*10 ceil(peaksM(2)/10)*10]);
        axis(ax(2), [0 triallen 0 100]);
        %set(ax(1), 'YTick', floor(peaksM(1)/10)*10:10:ceil(peaksM(2)/10)*10, 'FontSize', 8);
        %set(ax(1), 'YTickLabel', floor(peaksM(1)/10)*10:10:ceil(peaksM(2)/10)*10, 'FontSize', 8);
        set(ax(2), 'YTick', 0:20:40);
        set(ax(2), 'YTickLabel', 0:20:40);
    elseif (r{k}.chantype == 'O')  % all EMG channels have the same scaling
        axis(ax(1), [0 triallen floor(peaksO(1)/10)*10 ceil(peaksO(2)/10)*10]);
        axis(ax(2), [0 triallen 0 100]);
        %set(ax(1), 'YTick', floor(peaksO(1)/10)*10:10:ceil(peaksO(2)/10)*10, 'FontSize', 8);
        %set(ax(1), 'YTickLabel', floor(peaksO(1)/10)*10:10:ceil(peaksO(2)/10)*10, 'FontSize', 8);
        set(ax(2), 'YTick', 0:20:40);
        set(ax(2), 'YTickLabel', 0:20:40);
    else
        v = axis(ax(1));
        axis(ax(1), [0 triallen v(3) v(4)]);
        v = axis(ax(2));
        axis(ax(2), [0 triallen v(3) v(4)]);
    end;

    set(ax(1), 'XTick', 0:fs:triallen);
    set(ax(1), 'XTickLabel', startend(1):startend(2));
    set(ax(2), 'XTick', []);
    set(ax(2), 'XTickLabel', []);

    if (~isempty(cue) && cue >= startend(1))
        v = axis;
        line([(cue-startend(1))*fs (cue-startend(1))*fs], [v(3) v(4)], ...
             'Color', 'k');
    end;
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