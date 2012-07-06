% pop_ERPLAB_polydetrend() - interactively detrend EEG dataset data using ERPLAB_polydetrend()
%
% Usage:
%
%   >> [EEG, com] = pop_ERPLAB_polydetrend( EEG, window)
%
% Graphical interface:
%
%   "Window" - [edit box] minimun time window (in sec) to represent the
%                         mean of DC behaviour (into the window).
%                         Row data will divided per this window, generating
%                         a joint of points over which a polynomial will be
%                         calculated.
%                         Applying the same vector of time to polynomial we
%                         will get the full DC behavior per channel. The
%                         last one is subtracted from each channel finally.
%
% Inputs:
%
%   EEG       - input dataset
%   window    - minimun time window (in sec)
%
%
% Outputs:
%
%   EEGOUT    - output dataset
%
%     Author: ERPLAB Team, Center for Mind & Brain
%             Universidad de California, Davis. 2007


function [EEG, com] = pop_ERPLAB_polydetrend( EEG, window)

com = '';

if exist('filtfilt','file') ~= 2
    error('Warning: cannot find the signal processing toolbox')
end

if nargin < 1
    help pop_ERPLAB_polydetrend
    return
end

if isempty(EEG.data)
    disp('Pop_ERPLAB_polydetrend() error: cannot filter an empty dataset')
    return
end


uilist = { ...
    { 'style' 'text' 'string' ' Window length (s) - see help' } ...
    { 'style' 'edit' 'string' '' }};

geometry = { [1 0.35] };

result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Detrend the data -- pop_ERPLAB_polydetrend()', ...
    'helpcom', 'pophelp(''pop_ERPLAB_polydetrend'')');

if isempty(result) || isempty(result{1})
    return
end

mesg='Ups!!!';
orderpoly = 5;                          % polynomial order n : p1*x^n + p2*x^(n-1)+...+pn*x^1
windowsec = eval( result{1} );          % in sec
window = windowsec * (EEG.srate);       % sec converted to samples
numwin=round(EEG.pnts/window);          % number of points per window
durt= round(EEG.pnts/EEG.srate) / 60;   % recording duration en minutes
fprintf('\r Using %g windows of %g sec for a recording of %g min aprox.\n', numwin, windowsec, durt)


if numwin < orderpoly
    fprintf('***ERPLAB polydetrend Warning***  for %s', EEG.setname)
    fprintf(...
    '\r The number of estimated windows ( numwin = %g ) is lesser than the order of yhe polynomial (n = %g)our data.\n',...
    numwin, orderpoly)
    disp('decrease window please...')
    return
end

if windowsec < 3
    disp('Sorry. Your window is too close to ERP window...')
    disp('increase window please...')
    return
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_ERPLAB_polydetrend', EEG, 'warning', 'on', 'params', ...
        { window } );
    return
end

% Pensando en EEGOUT = ERPLAB_butter1( EEG, locutoff, hicutoff, order)
%[smoothdata,filtwts] = ERPLAB_butter1(data,srate,locutoff,hicutoff,epochframes,filterorder)
%options = { window };

%EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);

EEG.data = ERPLAB_polydetrend( EEG.data, window, orderpoly);

EEG.icaact = [];

%function [EEG, com] = pop_ERPLAB_butter1( EEG, window)
com = sprintf( '%s = pop_ERPLAB_polydetrend( %s, %s, %s);', inputname(1), inputname(1), ...
    num2str( window), num2str( orderpoly ));
return


