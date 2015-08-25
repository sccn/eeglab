function checkfunctionmatlab(func, toolbox)
if ismatlab
    octavefolder = fullfile(fileparts(which('eeglab.m')),'functions','octavefunc');
    if license('test',toolbox)
        pathtemp = which(func,'-all'); 
        if length(pathtemp) > 1
            fprintf( [ 'Octave functions should not run on Matlab' 10 'Removing path ... ' 10])
            hit = zeros(1,length(pathtemp));
            for i = 1 : length(pathtemp)
                pathtemp{i} = fileparts(pathtemp{i});
                tmp         = strfind(pathtemp{i},octavefolder);
                if ~isempty(tmp)
                    hit(i) = 1;
                end
            end
            cellfun(@(x) rmpath(x),pathtemp(find(hit)));   
        end
    else
        addpath(genpath(octavefolder));
        warning([ toolbox ' toolbox is absent or not in the path, using replacement functions' ]);
    end;
end;
