function plotSpectrum(r)

plot_index = r{1}.plot_index;
n_rows = r{1}.n_rows;
n_cols = r{1}.n_cols;
fs = r{1}.fs;
name = r{1}.name;
trials = r{1}.trials;
class = r{1}.class;
triallen = r{1}.triallen;
len_ref_int = r{1}.len_ref_int;
len_act_int = r{1}.len_act_int;
st_ref = r{1}.st_ref;
st_act = r{1}.st_act;
f = r{1}.f;

% Plot spectra
fig = figure('Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);
set(gcf, 'DefaultAxesFontSize', 6);

peaks = [0 0];  % Minimum and maximum

el_u = round((512 * f)/fs) + 1;  % Last frequency bin to plot

for k = 1:length(plot_index)
    subplot(n_rows, n_cols, plot_index(k));

    sp_ref = [];
    for l = 1:trials
        temp = fft(st_ref(k, (l-1)*len_ref_int + 1:l*len_ref_int).*...
                   hanning(len_ref_int)', 512);
        sp_ref = [sp_ref (temp.*conj(temp)/512)'];
    end;
    sp_ref = mean(sp_ref,2);
    semilogy(fs*(0:el_u-1)/512, sp_ref(1:el_u));

    hold on;

    sp_act = [];
    for l = 1:trials
        temp = fft(st_act(k, (l-1)*len_act_int + 1:l*len_act_int).*hanning(len_act_int)', 512);
        sp_act = [sp_act (temp.*conj(temp)/512)'];
    end;
    sp_act = mean(sp_act,2);
    semilogy(fs*(0:el_u-1)/512, sp_act(1:el_u), 'r');

    if (min(sp_ref) < peaks(1))
        peaks(1) = min(sp_ref);
    end;
    if (max(sp_ref) > peaks(2))
        peaks(2) = max(sp_ref);
    end;
    if (min(sp_act) < peaks(1))
        peaks(1) = min(sp_act);
    end;
    if (max(sp_act) > peaks(2))
        peaks(2) = max(sp_act);
    end;

    set(gca, 'XTick', 0:20:f);
    set(gca, 'XTickLabel', 0:20:f);
    set(gca, 'YTickLabel', []);
end;

peaks(1) = floor(peaks(1)/10)*10;
peaks(2) = ceil(peaks(2)/10)*10;

for k = 1:length(plot_index)
    subplot(n_rows, n_cols, plot_index(k));
    axis([0 f peaks]);
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