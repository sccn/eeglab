% iseeglabdeployed - true for EEGLAB compile version and false otherwise

function val = iseeglabdeployed;
%val = 1; return;
if exist('isdeployed')
     val = isdeployed;
else val = 0;
end;
