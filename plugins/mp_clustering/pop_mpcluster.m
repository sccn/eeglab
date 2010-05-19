% pop_mpcluster() - Clusters STUDY ICs using Measure Product method. 
%                   In this method, IC measures, except equiv. dipoles, (ERP, ERSP...)  
%                   are compared for each IC pair and their dissimilarity is multiplied
%                   together to form a combined pairwise dissimilarity matrix. This matrix
%                   is then normalized, weighted and added to the normalized and weighted 
%                   IC equiv. dipole distance matrix. The final dissimilarity matrix is
%                   then clustered using affinity clustering  method. 
%                   You can control the effect of equiv. dipole distances in
%                   the clustering by setting the 'Relative dipole weight'
%                   parameter in the pop-uo GUI. For example, by setting
%                   this value to 0.8, the final dissimilarity matrix will consist of 80% 
%                   distance dissimilarity and 20% of other measures combined together.
%                   Please note that the number of returned clusters may slighty 
%                   differ from the number requested in the GUI.
%                  
%
% Usage:
%     >> STUDY = pop_mpcluster(STUDY, ALLEEG) % popup window
%
% See also:  std_mpcluster(), std_mpreclust(), , pop_mpreclust(), apclusterK()
% 
% Author: Nima Bigdely-Shamlo, SCCN/INC/UCSD, 2009

function [STUDY ALLEEG command] = pop_mpcluster(STUDY, ALLEEG)
% command is used for keeping a history.
% disable measure checkboxes which are not present (calculated) in
% pre-clustering.

% ERP
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'erpCorr')
    erpEnable = 'on';
else
    erpEnable = 'off';
end;
erpChecked = strcmp(erpEnable, 'on'); % only check items if they are enabled.

% ERSP
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'erspCorr')
    erspEnable = 'on';
else
    erspEnable = 'off';
end;
erspChecked = strcmp(erspEnable, 'on'); % only check items if they are enabled.

% ITC
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'itcCorr')
    itcEnable = 'on';
else
    itcEnable = 'off';
end;
itcChecked = strcmp(itcEnable, 'on'); % only check items if they are enabled.

% dipole
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'compDistance')
    dipoleEnable = 'on';
else
    dipoleEnable = 'off';
end;
dipoleChecked = strcmp(dipoleEnable, 'on'); % only check items if they are enabled.


% spectra
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'specCorr')
    specEnable = 'on';
else
    specEnable = 'off';
end;
specChecked = strcmp(specEnable, 'on'); % only check items if they are enabled.


% scalp map
if isfield(STUDY.preclust,'similarity') && isfield(STUDY.preclust.similarity,'mapCorr')
    scalpEnable = 'on';
else
    scalpEnable = 'off';
end;

% popup GUI

returnedFromGui = inputgui( 'geometry', { [1 0.5]  [1 0.5] 1 1 [1 1 1] [1 1 1] [1 1 1] [1 1 1] [1 1 1] [1 1 1] 1 [3 1.5]}, ...
    'geomvert', [], 'uilist', { ...
    { 'style', 'text', 'string', 'Number of clusters to compute:'}, ...
    { 'style', 'edit', 'string', '10' 'tag' 'numberOfClusters' } , ...
       { 'style', 'text', 'string', 'Relative dipole weight (between 0 and 1):'}, ...
    { 'style', 'edit', 'string', '0.8' 'tag' 'numberOfClusters' } , ...
        { 'style', 'text', 'string', [ 'Select measuretures to be used in the clustering:' ] }, {}, ...
    {},{ 'Style', 'checkbox', 'string' 'Dipole' 'tag' 'scale' 'value' 1} , {}, ...
    {},{ 'Style', 'checkbox', 'string' 'ERP' 'tag' 'scale' 'value' erpChecked 'enable' erpEnable}, {},...
    {},{ 'Style', 'checkbox', 'string' 'ERSP' 'tag' 'scale' 'value' erspChecked 'enable' erspEnable}, {},...
    {},{ 'Style', 'checkbox', 'string' 'ITC' 'tag' 'scale' 'value' itcChecked 'enable' itcEnable}, {},...
    {},{ 'Style', 'checkbox', 'string' 'Spectra' 'tag' 'scale' 'value' specChecked 'enable' specEnable}, {}, ...
    {}, { 'Style', 'checkbox', 'string' 'Scalp map' 'tag' 'scale' 'value' 0 'enable' scalpEnable}, {},...
    {}, { 'Style', 'checkbox', 'string' 'Separate outliers (enter std.)' 'tag' 'scale' 'value' 1}, { 'style', 'edit', 'string', '3' 'tag' 'outlierSTD' }, ...

    }, 'helpcom','pophelp(''pop_mpcluster'');', 'title', 'Measure Product clustering -- pop_mpcluster()');



if isempty(returnedFromGui) % an empty returnedFromGui means the Cancel button has been pressed so nothing should be done.
    command = '';
    return; % Cancel button is pressed, so do nothing.
else
    
    % analysze answers returned from the GUI
    numberOfClusters = str2num(returnedFromGui{1});
    methodParameter = str2num(returnedFromGui{2});
    
    answers = cell2mat(returnedFromGui(3:end-1));
    measureNamesInGUIorder = {'dipole', 'erp', 'ersp', 'itc', 'spec', 'map'};
    measuresToUseInClustering = measureNamesInGUIorder(find(answers(1:end-1))); %#ok<FNDSB>
    
    if answers(end) % checkbox for outlier
        outlierSTD = str2num(returnedFromGui{end});
    else
        outlierSTD = Inf;
    end;
    
    STUDY = std_mpcluster(STUDY, ALLEEG, numberOfClusters, outlierSTD, measuresToUseInClustering, methodParameter);
    
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
    
    % pop up the cluster edit and visualization .
    [STUDY commandFromPop_clustedit] = pop_clustedit(STUDY, ALLEEG); 
    
    command = ['STUDY = std_mpcluster(STUDY, ALLEEG, ' num2str(numberOfClusters) ', ' num2str(outlierSTD) ', {' measuresInOneString '} , ' num2str(methodParameter) ');'];
    command = [command '\n' commandFromPop_clustedit];     % add the command from pop_clustedit() to the history too.
    STUDY.history =  sprintf('%s\n%s',  STUDY.history, command);
        
end;