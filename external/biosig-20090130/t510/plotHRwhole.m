function plotHRwhole(r)

hr = r{1}.hr;
fs = r{1}.fs;
startend = r{1}.startend;
h = r{1}.h;

fig = figure('Color', [1 1 1]);

t = 0:1/fs:length(hr)/fs-1/fs;
plot(t, hr);
hold on;

if (isempty(startend))
    startend = [0 t(end)];
end;

% Determine minimum and maximum of y-axis
min_hr = floor(min(hr) * 10) / 10;
max_hr = ceil(max(hr) * 10) / 10;

axis([startend(1) startend(2) min_hr max_hr]);

% Mark cues as vertical lines in different colors
% At the moment, the following event types are recognized as cues:
% 0x301, 0x302, 0x303, 0x304, 0x306

cue = h.EVENT.POS(h.EVENT.TYP==hex2dec('301') & h.EVENT.POS>=startend(1)*fs & h.EVENT.POS<=startend(2)*fs)/fs;
if (~isempty(cue))
    line([cue cue], [min_hr max_hr], 'Color', 'k');
end;

cue = h.EVENT.POS(h.EVENT.TYP==hex2dec('302') & h.EVENT.POS>=startend(1)*fs & h.EVENT.POS<=startend(2)*fs)/fs;
if (~isempty(cue))
    line([cue cue], [min_hr max_hr], 'Color', 'r');
end;

cue = h.EVENT.POS(h.EVENT.TYP==hex2dec('303') & h.EVENT.POS>=startend(1)*fs & h.EVENT.POS<=startend(2)*fs)/fs;
if (~isempty(cue))
    line([cue cue], [min_hr max_hr], 'Color', 'g');
end;

cue = h.EVENT.POS(h.EVENT.TYP==hex2dec('304') & h.EVENT.POS>=startend(1)*fs & h.EVENT.POS<=startend(2)*fs)/fs;
if (~isempty(cue))
    line([cue cue], [min_hr max_hr], 'Color', 'c');
end;

cue = h.EVENT.POS(h.EVENT.TYP==hex2dec('306') & h.EVENT.POS>=startend(1)*fs & h.EVENT.POS<=startend(2)*fs)/fs;
if (~isempty(cue))
    line([cue cue], [min_hr max_hr], 'Color', 'm');
end;

% Set suitable paper size
set(gcf, 'Color', [1 1 1]);
set(gcf, 'PaperOrientation', 'landscape');
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperType', 'A4');
set(gcf, 'PaperPosition', [1, 1, 27.7, 19]);