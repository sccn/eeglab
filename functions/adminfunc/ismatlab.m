function res = ismatlab
% true if called from Matlab; false if called from Octave

v = version;
if v(1) > '6'
    res = 1;
else
    res = 0;
end
