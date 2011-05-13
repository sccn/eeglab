function res = ismatlab;

v = version;
if v(1) > '4'
    res = 1;
else
    res = 0;
end;