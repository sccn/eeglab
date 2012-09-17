% subfunction checking the number of local copies
% -----------------------------------------------
function ncopies = checkcopies_local(obj, arg);
ncopies = 0;
if isstruct(arg)
    for ilen = 1:length(arg)
        for index = fieldnames(arg)'
            ncopies = ncopies + checkcopies_local(obj, arg(ilen).(index{1}));
            if ncopies > 1, return; end;
        end;
    end;
elseif iscell(arg)
    for index = 1:length(arg(:))
        ncopies = ncopies + checkcopies_local(obj, arg{index});
        if ncopies > 1, return; end;
    end;
elseif isa(arg, 'mmo') && isequal(obj, arg)
    ncopies = 1;
else
    ncopies = 0;
end;