% pop_ERPLAB_butter1() - interactively filter EEG dataset data using ERPLAB_butter1()
%
% Usage:
%
%   >> [EEG, com] = pop_ERPLAB_butter1( EEG, locutoff, hicutoff, filterorder)
%
% Graphical interface:
%   "Lower edge ..." - [edit box] Lower edge of the frequency pass band (Hz)
%                 Same as the 'locutoff' command line input.
%   "Higher edge ..." - [edit box] Higher edge of the frequency pass band (Hz)
%                 Same as the 'hicutoff' command line input.
%
% Inputs:
%
%   EEG       - input dataset
%   locutoff  - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%   hicutoff  - higher edge of the frequency pass band (Hz) {0 -> highpass}
%   filterorder - length of the filter in points {default 3*fix(srate/locutoff)}
%
%
% Outputs:
%
%   EEGOUT    - output dataset
%
%     Author: ERPLAB Team, Center for Mind & Brain
%             Universidad de California, Davis. 2007
%
%  Modified error with inputgui for multiple file filtering. JLC Dic 7,2007


function [EEG, com] = pop_ERPLAB_butter1( EEG, locutoff, hicutoff, filterorder)

com = '';

if exist('filtfilt','file') ~= 2
    error('Warning: cannot find the signal processing toolbox')
end

if nargin < 1
    help pop_ERPLAB_butter1
    return
end

if isempty(EEG(1).data)
    disp('Pop_ERPLAB_butter1() error: cannot filter an empty dataset')
    return
end


if nargin < 2   %*******  Dic 7, 2007

    uilist = { ...
        { 'style' 'text' 'string' 'High Pass Cutoff (Hz)' } ...
        { 'style' 'edit' 'string' '0.05' } ...
        { 'style' 'text' 'string' 'FIR Filter order (default is 5)' } ...
        { 'style' 'edit' 'string' '' }};

    geometry = { [3 1] [3 1] };

    result = inputgui( 'geometry', geometry, 'uilist', uilist, 'title', 'Filter the data -- pop_ERPLAB_butter1()', ...
        'helpcom', 'pophelp(''pop_ERPLAB_butter1'')');

    if isempty(result)
        return
    end

    if isempty(result{1})
        result{1} = '0';
    end

    locutoff = eval( result{1} );
    hicutoff = 0;

    if isempty( result{2} )
        filterorder = 5;
    else
        filterorder = eval( result{2} );
    end


else  %*******  Dic 7, 2007
    if nargin < 3
        hicutoff = 0;
    end
    if nargin < 4
        filtorder = 5;
    end
   
end


if locutoff == 0 && hicutoff == 0
    disp('what happened with you?')
    return
end

% process multiple datasets
% -------------------------
if length(EEG) > 1
    [ EEG com ] = eeg_eval( 'pop_ERPLAB_butter1', EEG, 'warning', 'on', 'params', ...
        { locutoff, hicutoff, filterorder} );
    return
end

% Pensando en EEGOUT = ERPLAB_butter1( EEG, locutoff, hicutoff, order)
%[smoothdata,filtwts] = ERPLAB_butter1(data,srate,locutoff,hicutoff,epochframes,filterorder)
options = { EEG.srate, locutoff, hicutoff, EEG.pnts, filterorder};

%EEG.data = reshape(EEG.data, EEG.nbchan, EEG.pnts*EEG.trials);

EEG.data = ERPLAB_butter1( EEG.data, options{:});

EEG.icaact = [];

com = sprintf( '%s = pop_ERPLAB_butter1( %s, %s, %s, %s );', inputname(1), inputname(1), ...
    num2str( locutoff), num2str( hicutoff), num2str( filterorder ));
return


