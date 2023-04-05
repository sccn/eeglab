function nameout = cleanvarname(namein)

% Routine to remove unallowed characters from strings
% nameout can be use as a variable or field in a structure
% Cyril Pernet & Arnaud Delorme

% custom change
if strcmp(namein,'#')
    namein = 'nb';
end

% 1st change space to underscore (keeps readability)
for l=1:length(namein)
    if isspace(namein(l))
        namein(l) = '_';
    end
end

% keep letters, numbers and underscore
keep        = isstrprop(namein, 'alpha') + isstrprop(namein, 'digit');
underscores = zeros(1,size(namein,2));
underscores(strfind(namein,'_')) = 1;
keep    = keep + underscores;
nameout = namein(logical(keep));

% remove usual suspects JIC isstrprop did not work
stringcheck = [strfind(nameout,'('), ...
    strfind(nameout,')') ...
    strfind(nameout,'[') ...
    strfind(nameout,']') ...
    strfind(nameout,'{') ...
    strfind(nameout,'}') ...
    strfind(nameout,'-') ...
    strfind(nameout,'+') ...
    strfind(nameout,'*') ...
    strfind(nameout,'/') ...
    strfind(nameout,'#') ...
    strfind(nameout,'%') ...
    strfind(nameout,'&') ...
    strfind(nameout,'@') ...
    ];

if ~isempty(stringcheck)
    nameout(stringcheck) = [];
end

% No number in position 1
nb_check = 1;
while nb_check
    if isletter(nameout(1)) 
       nb_check = 0;
    elseif isnumeric(eval(nameout(1)))
        nameout(1) = [];
    end
end

% last check
if ~isvarname(nameout)
    error('the variable name to use is still invalid, check chars to remove')
end
