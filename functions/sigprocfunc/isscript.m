% function checking if a specific file is a script%

function bool = isscript(fileName);

fid = fopen(fileName, 'r');
cont = true;
while cont && ~feof(fid)
    l = strtok(fgetl(fid));
    
    if ~isempty(l) && l(1) ~= '%'
        if strcmpi(l, 'function')
            bool = false; return;
        else
            bool = true; return;
        end;
    end;
end;
bool = true; return;
