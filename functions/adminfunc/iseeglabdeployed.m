% iseeglabdeployed - true for EEGLAB compile version and false otherwise
function val = iseeglabdeployed
try
     val = isdeployed;
catch
    val = 0;
end