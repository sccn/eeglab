
% eeg_toolbox_doc_make - script to use m2html to create eeg_toolbox documentation

if exist('m2html'),
    location = which('eeg_toolbox');
    [path,file,ext] = fileparts(location);
    path = strrep(path,'eeg_toolbox','');
    cd(path);
    m2html('m','eeg_toolbox','html',['eeg_toolbox',filesep,'doc_m2html']);
else
    error('cannot locate m2html, see http://www.madic.org/download/matlab/m2html/');
end

clear ext file location path
