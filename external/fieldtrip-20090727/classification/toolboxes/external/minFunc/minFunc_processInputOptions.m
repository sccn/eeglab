
function [verbose,verboseI,debug,doPlot,maxFunEvals,maxIter,tolFun,tolX,method,...
    corrections,c1,c2,LS_init,LS,cgSolve,SR1,cgUpdate,initialHessType,...
    HessianModify,Fref,useComplex,numDiff,LS_saveHessianComp,...
    DerivativeCheck,Damped,HvFunc] = ...
    minFunc_processInputOptions(o)
verbose = 1;
verboseI= 1;
debug = 0;
doPlot = 0;
method = 3;
cgSolve = 0;

o = toUpper(o);

if isfield(o,'DISPLAY')
    switch(upper(o.DISPLAY))
        case 0
            verbose = 0;
            verboseI = 0;
        case 'FINAL'
            verboseI = 0;
        case 'OFF'
            verbose = 0;
            verboseI = 0;
        case 'NONE'
            verbose = 0;
            verboseI = 0;
        case 'FULL'
            debug = 1;
        case 'EXCESSIVE'
            debug = 1;
            doPlot = 1;
    end
end


if isfield(o,'METHOD')
    m = upper(o.Method);
    switch(m)
        case 'TENSOR'
            method = 8;
        case 'NEWTON'
            method = 7;
        case 'NEWTON0LBFGS'
            method = 5;
            cgSolve = 1;
        case 'NEWTON0'
            method = 5;
        case 'BFGS'
            method = 4;
        case 'LBFGS'
            method = 3;
        case 'BB'
            method = 2;
        case 'CG'
            method = 1;
        case 'SD'
            method = 0;
    end
end

c2 = 0.9;
LS_init = 0;
LS = 4;
Fref = 1;
% Method Specific Default Options if different than the above
if method == 2
    LS = 2;
    Fref = 10;
elseif method == 1
    c2 = 0.2;
    LS_init = 1;
elseif method == 0
    LS = 2;
    LS_init = 1;
end

maxFunEvals = getOpt(o,'MAXFUNEVALS',1000);
maxIter = getOpt(o,'MAXITER',500);
tolFun = getOpt(o,'TOLFUN',1e-5);
tolX = getOpt(o,'TOLX',1e-9);
corrections = getOpt(o,'CORR',100);
c1 = getOpt(o,'C1',1e-4);
c2 = getOpt(o,'C2',c2);
LS_init = getOpt(o,'LS_INIT',LS_init);
LS = getOpt(o,'LS',LS);
cgSolve = getOpt(o,'CGSOLVE',cgSolve);
SR1 = getOpt(o,'SR1',0);
cgUpdate = getOpt(o,'CGUPDATE',1);
initialHessType = getOpt(o,'INITIALHESSTYPE',1);
HessianModify = getOpt(o,'HESSIANMODIFY',0);
Fref = getOpt(o,'FREF',Fref);
useComplex = getOpt(o,'USECOMPLEX',0);
numDiff = getOpt(o,'NUMDIFF',0);
LS_saveHessianComp = getOpt(o,'LS_SAVEHESSIANCOMP',1);
DerivativeCheck = getOpt(o,'DERIVATIVECHECK',0);
Damped = getOpt(o,'DAMPED',1);
HvFunc = getOpt(o,'HVFUNC',[]);
end

function [v] = getOpt(options,opt,default)
if isfield(options,opt)
    if ~isempty(getfield(options,opt))
        v = getfield(options,opt);
    else
        v = default;
    end
else
    v = default;
end
end

function [o] = toUpper(o)
if ~isempty(o)
    fn = fieldnames(o);
    for i = 1:length(fn)
        o = setfield(o,upper(fn{i}),getfield(o,fn{i}));
    end
end
end