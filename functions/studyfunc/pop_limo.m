function [STUDY,com] = pop_limo(STUDY, ALLEEG, varargin)

if nargin < 2
    help pop_limo;
    return;
end;
com = '';
if nargin < 3
    dataMeasures = { 'ERP' 'Spectrum' };
    fileMeasures = { 'daterp' 'datspec' };
    methods      = { 'OLS' 'WLS' };

    uilist = { ...
        {'style' 'text'       'string' [ 'LInear MOdeling of EEG data of design' int2str(STUDY.currentdesign) ] 'fontweight' 'bold' 'fontsize', 12} ...
        {'style' 'text'       'string' [ '(Use STUDY design interface to switch to a different design)' ] } {} ...
        {'style' 'text'       'string' 'Input data to use for the GLM' } ...
        {'style' 'popupmenu'  'string' dataMeasures 'tag' 'measure' } ...
        {'style' 'text'       'string' 'Optimization method' } ...
        {'style' 'popupmenu'  'string' methods 'tag' 'method' } ...
        {'style' 'checkbox'   'string' 'Erase all previous models for this design' 'tag' 'erase' } };
    
    cline = [1.1 1.1];
    geometry = { 1 1 1 cline cline 1 };
    geomvert = [ 1 1 1 1     1     1 ];
        
    [out_param userdat tmp res] = inputgui( 'geometry' , geometry, 'uilist', uilist, 'geomvert', geomvert, ...
                                            'title', 'LInear MOdeling of EEG data -- pop_limo()');
    if isempty(res), return; end;
    
    options = { 'method' methods{res.method} 'measure' fileMeasures{res.measure} 'erase' fastif(res.erase, 'on', 'off') };
else
    options = varargin;
end;
if strcmp(fastif(res.measure, 'on', 'off'), 'on')
    if exist([STUDY.filepath filesep 'limo_batch_report'],'dir') == 7
        rmdir([STUDY.filepath filesep 'limo_batch_report'],'s');
    end
end
[tmp,STUDY] = std_eeglab2limo(STUDY, ALLEEG, options{:});
com = sprintf('std_eeglab2limo(STUDY, ALLEEG, %s);', vararg2str(options));
