function val = iseeglabdeployed
% iseeglabdeployed - true for EEGLAB compile version and false otherwise
try
     val = isdeployed;
catch
    val = 0;
end
