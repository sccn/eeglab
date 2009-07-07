function [x,f,exitflag,output] = minFunc(funObj,x0,Aeq,Beq,options,varargin)
% minFunc(funObj,x0,Aeq,Beq,options,varargin)
%
% Wrapper that lets you use minFunc to minimize subject to 
%   linear equality constraints: Aeq*x == Beq
% Aeq(constraint,variable)
% beq(constraint,1)
%
% Todo: Take advantage of sparsity

[verbose,verboseI,debug,doPlot,maxFunEvals,maxIter,tolFun,tolX,method,...
    corrections,c1,c2,LS_init,LS,cgSolve,SR1,cgUpdate,initialHessType,...
    HessianModify,Fref,useComplex,numDiff,LS_saveHessianComp,DerivativeCheck] = minFunc_processInputOptions(options);


% Check that A has the right size
p = length(x0);
if p ~= size(Aeq,2)
    fprintf('Number of columns in Aeq not consistent with x0\n');
    pause;
    return;
end

% Compute null-space of A by orthogonal decomposition
m = size(Aeq,1);
[Q,R] = qr(Aeq');
Z = Q(1:p,(m+1):p);

% Check that starting point is feasible
if Aeq*x0 == Beq
    if debug
        fprintf('Starting Point is Feasible\n');
    end
else
    % If not, compute a feasible starting point
    if debug
        fprintf('Starting Point Not Feasible\n');
    end
    AQ = Aeq*Q;
    x0 = Q(1:p,1:m)*(AQ(1:m,1:m)\Beq);
end

varargin = {Z,funObj,varargin{:}};
funObj = @ReducedFunObj;

% Only backtracking line-searches can be used
options.LS = min(LS,2);

% Solve
[x,f,exitflag,output] = minFunc(funObj,x0,options,varargin{:});
end

function [f,g,H] = ReducedFunObj(x,Z,funObj,varargin)

if nargout == 1
    f = funObj(x,varargin{:});
elseif nargout == 2
    [f,g] = funObj(x,varargin{:});
    g = Z'*g;
else
    [f,g,H] = funObj(x,varargin{:});
    g = Z'*g;
    H = Z'*H*Z;
end

end