function ver = eegplugin_BCI2000import(fig, trystrs, catchstrs)

ver = 'BCI2000import0.2';
if nargin < 3
    error('eegplugin_loadBCI2000 requires 3 arguments');
end;

% find import data menu
% ---------------------
menu = findobj(fig, 'tag', 'import data');

% menu callbacks
% --------------
comcnt = [ trystrs.no_check 'EEG = pop_loadBCI2000();'     catchstrs.new_and_hist ];

% create menus
% ------------
uimenu( menu, 'label', 'From BCI2000 .DAT file', 'callback', comcnt, 'separator', 'on' );
