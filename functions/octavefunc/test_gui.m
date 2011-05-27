% This simple scripts test Octave GUI capabilities
% (it crashes Octave under OS X)

% The line is sometimes necessary
graphics_toolkit fltk 

figure;
tmp  = uimenu('label', 'This is a custom menu', 'callback', 'disp(''Hello'');');
tmp2 = uimenu(tmp, 'label', 'Menu showing "Hello"', 'callback', 'disp(''Hello'');');
tmp2 = uimenu(tmp, 'label', 'Menu plotting a new figure', 'callback', 'figure; plot([1:10]);');
