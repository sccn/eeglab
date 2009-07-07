function [x,f,exitflag,output] = minFunc(funObj,x0,A,B,options,varargin)
% minFunc(funObj,x0,A,B,options,varargin)
%
% Linear inequality-constrained optimizer using a Log-Barrier
%
% minimizes funObj subject to A*x <= B
%
% set options.cont = 0 to use a 'short-cut' strategy (mu fixed)
% set options.cont = 1 to use a 'continuation' strategy (mu increased)
%
% x0 must be feasible
%
% Todo: Implement a better way to handle inequality constraints

% Get Parameters
[verbose,verboseI,debug,doPlot,maxFunEvals,maxIter,tolFun,tolX,method,...
    corrections,c1,c2,LS_init,LS,cgSolve,SR1,cgUpdate,initialHessType,...
    HessianModify,Fref,useComplex,numDiff,LS_saveHessianComp,DerivativeCheck] = minFunc_processInputOptions(options);

% Initialize
p = length(x0);
x = x0;
d = zeros(p,1);
t = 1;

% Test Presence of Mex Files
if (exist('lbfgsC','file')==3)
    lbfgsDir = @lbfgsC;
else
    lbfgsDir = @lbfgs;
end

% Form Barrier-Augmented Objective Function
varargin = {A,B,funObj,varargin{:}};
funObj = @barrierFunObj;

% Only backtracking line-searches can be used
options.LS = min(LS,2);

if isfield(options,'mu') && ~isempty(options.mu)
    mu = options.mu;
else
    mu = 1e-7;
end

if isfield(options,'cont') && options.cont==1
    muCurr = 1;
    x = x0;
    while muCurr > mu
        [x,f,exitflag,output] = minFunc(funObj,x,options,muCurr,varargin{:});
        muCurr = muCurr/10;
    end
else
    [x,f,exitflag,output] = minFunc(funObj,x0,options,mu,varargin{:});
end

end

function [f,g,H] = barrierFunObj(x,mu,A,B,funObj,varargin)
[m,p] = size(A);

if nargout == 1
    f = funObj(x,varargin{:});
    f = f - mu*sum(log(B-A*x));
elseif nargout == 2
    [f,g] = funObj(x,varargin{:});
    f = f - mu*sum(log(B-A*x));
    g = g + mu*sum(A.'./repmat((B-A*x).',[p 1]),2);
else
    [f,g,H] = funObj(x,varargin{:});
    f = f - mu*sum(log(B-A*x));
    g = g + mu*sum(A.'./repmat((B-A*x).',[p 1]),2);
    H = H + mu*A.'*(A./repmat((B-A*x).^2,[1 p]));
end

end