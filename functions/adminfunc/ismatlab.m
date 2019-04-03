% true if called from Matlab; false if called from Octave

function res = ismatlab

v = version;
if v(1) > '6'
    res = 1;
else
    res = 0;
end
