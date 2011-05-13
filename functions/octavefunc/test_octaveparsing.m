function test_octaveparsing;

a = dir;

for index = 1:length(a)
    %fprintf('testing %s\n', a(index).name);
    if exist(a(index).name) == 7 && ~strcmpi(a(index).name, '.') && ~strcmpi(a(index).name, '..')
        cd(a(index).name);
        test_octaveparsing;
        cd('..');
    elseif length(a(index).name) > 2 && a(index).name(end) == 'm'
        try
            which(a(index).name(1:end-2));
        catch
            tmp = lasterr;
            inds = find(tmp == 10);
            fprintf('*************\n%s\n*************\n', tmp(1:inds(1)-1));
        end;
    end;
end;