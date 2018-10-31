function electrodes_to_tsv(EEG)

% From an EEG variable (i.e. EEG=pop_loadset(*.set), export the EEG.channel
% locations as tsv file and initiate the json file following the BIDS 
% specification
%
% FORMAT electrodes_to_tsv(EEG)
%
% Author: Cyril Pernet - LIMO Team, University of Edinurgh

% electrode.tsv
for electrode = 1:size(EEG.chanlocs,2)
    ename{electrode} = EEG.chanlocs(electrode).labels; 
    x(electrode)     = EEG.chanlocs(electrode).X;
    y(electrode)     = EEG.chanlocs(electrode).Y;
    z(electrode)     = EEG.chanlocs(electrode).Z;
    type{electrode}  = EEG.chanlocs(electrode).type;
end

t = table(ename',x',y',z',type','VariableNames',{'name','x','y','z','type'});
electrodes_tsv_name = [EEG.filepath filesep EEG.filename(1:strfind(EEG.filename,'run-')-1) 'electrodes.tsv'];
writetable(t,electrodes_tsv_name,'FileType','text','Delimiter','\t');

% coordsystem.json
EEGCoordinateSystem = 'RAS';
EEGCoordinateUnits  = 'mm';
json = struct('EEGCoordinateSystem','RAS', ...
        'EEGCoordinateUnits','mm');
jsonwrite([EEG.filepath filesep EEG.filename(1:strfind(EEG.filename,'run-')-1) 'electrodes.json'],json,struct('indent','  '))


