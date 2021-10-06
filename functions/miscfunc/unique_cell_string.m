function uniqueStrings = unique_cell_string(c)
% unique string from a cell-array containing only strings, ignores all
% non-strings.

nonStringCells = [];
for i=1:length(c) % remove non-string cells
    if ~strcmp(class(c{i}),'char')
        nonStringCells = [nonStringCells i];
    end
end
c(nonStringCells) = [];

uniqueStrings = {};
for i=1:length(c) % remove non-string cells
    if ~isAlreadyEncountered(c{i}, uniqueStrings);
        uniqueStrings{end+1} = c{i};
    end
end

function result = isAlreadyEncountered(s, u) 
result = false;
for i=1:length(u)
    if strcmp(u{i}, s)
        result = true;
    end
end
        
    
