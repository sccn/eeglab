% remove all path with a given parent path

% varargin contains a list of path to exclude
function removepath(parentpath, varargin)

if isempty(parentpath), return; end;
folder = path;
if ispc, sep = ';'; else sep = ':'; end;
indSep = find(folder == sep);

indSep = [ 0 indSep length(folder)+1 ];
for iSep = 1:length(indSep)-1
    curPath = folder(indSep(iSep)+1:indSep(iSep+1)-1);
    if ~isempty(strfind(curPath, parentpath)) && ~any(strmatch(curPath, varargin, 'exact'))
        rmpath(curPath);
    end;
end;
