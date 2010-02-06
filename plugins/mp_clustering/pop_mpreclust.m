% pop_mpreclust() - Calculates pairwise similarity matrices for select EEG
%                   measures (erp, ersp..) to be used (later) in MP clustering.
%                   Pre-clustering results are placed under STUDY.preclust.similarity 
%                   field. This functions calls std_mpreclust() internally.
%
% Usage:
%     >> STUDY = pop_mpreclust(STUDY, ALLEEG) % popup window
%
% See also:  std_mpreclust(), std_mpcluster(), pop_mpcluster()
% 
% Author: Nima Bigdely-Shamlo, SCCN/INC/UCSD, 2009

function [STUDY, ALLEEG, command] = pop_mpreclust(STUDY, ALLEEG)
% command is used for keeping a history.
returnedFromGui = inputgui( 'geometry', { 1 1 [1 2 1] [1 2 1] [1 2 1] [1 2 1] [1 2 1] [1  2  1] 1 1 }, ...
    'geomvert', [], 'uilist', { ...
    { 'style', 'text', 'string', [ 'Select measures for Measure Product pre-clustering:' ] }, {}, ...
    {},{ 'Style', 'checkbox', 'string' 'Equiv. dipoles' 'tag' 'scale' 'value' 1} , {}, ...
    {},{ 'Style', 'checkbox', 'string' 'ERPs' 'tag' 'scale' 'value' 1}, {},...
    {},{ 'Style', 'checkbox', 'string' 'ERSPs' 'tag' 'scale' 'value' 1}, {},...
    {},{ 'Style', 'checkbox', 'string' 'ITCs' 'tag' 'scale' 'value' 1}, {},...
    {},{ 'Style', 'checkbox', 'string' 'Mean spectra' 'tag' 'scale' 'value' 1}, {}, ...
    {}, { 'Style', 'checkbox', 'string' 'Scalp maps' 'tag' 'scale' 'value' 0}, {},...
    {}, { 'Style', 'checkbox', 'string' 'Re-calculate All' 'tag' 'scale' 'value' 0}, ...
    }, 'helpcom','pophelp(''pop_mpreclust'');', 'title', 'MP pre-clustering -- pop_mpreclust()');

if isempty(returnedFromGui) % an empty returnedFromGui means the Cancel button has been pressed so nothing should be done.
    command = '';
    return; % Cancel button is pressed, so do nothing.
else
    
    answers = cell2mat(returnedFromGui);
    measureNamesInGUIorder = {'dipole', 'erp', 'ersp', 'itc', 'spec', 'map'};
    
    measuresToUseInClustering = measureNamesInGUIorder(find(answers(1:end-1))); %#ok<FNDSB>
    reCalculateAll = answers(end);
    STUDY = std_mpreclust(STUDY,ALLEEG, measuresToUseInClustering, reCalculateAll);
    
    % prepare 'command' variable for placing both in eeglab histry (accessible with eegh() ) and also
    % adding to  STUDY.history
    
    measuresInOneString = [];
    for i=1:length(measuresToUseInClustering)
        if i>1
            measuresInOneString = [measuresInOneString ' , ' '''' measuresToUseInClustering{i} ''''];
        else
            measuresInOneString = ['''' measuresToUseInClustering{1} ''''];
        end;
    end;
    
    command = ['[STUDY ALLEEG]= std_mpreclust(STUDY,ALLEEG, {' measuresInOneString '} , ' num2str(reCalculateAll) ');'];
    STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);
end;