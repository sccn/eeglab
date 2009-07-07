function [x,f,exitflag,output] = minFunc(funObj,x0,LB,UB,options,varargin)
% minFunc(funObj,x0,LB,UB,options,varargin)
%
% Bound-constrained optimizer
%
% The input/output format is the same as minFunc, with the following
%   additional args:
% LB(var,1): lower-bound on each variable (set to -inf for unbounded)
% UB(var,1): upper-bound on each variable (set be inf for unbounded)
%
% Todo: Fix Armijo condition along projection arc

% Get Parameters
[verbose,verboseI,debug,doPlot,maxFunEvals,maxIter,tolFun,tolX,method,...
    corrections,c1,c2,LS_init,LS,cgSolve,SR1,cgUpdate,initialHessType,...
    HessianModify,Fref,useComplex,numDiff,LS_saveHessianComp,DerivativeCheck] = minFunc_processInputOptions(options);

% Initialize
p = length(x0);
d = zeros(p,1);
t = 1;

% Test Presence of Mex Files
if (exist('lbfgsC','file')==3)
    lbfgsDir = @lbfgsC;
else
    lbfgsDir = @lbfgs;
end

% If necessary, form numerical differentiation functions
if numDiff
    varargin(3:end+2) = varargin(1:end);
    varargin{1} = useComplex;
    varargin{2} = funObj;
    if method < 7
        if debug
            fprintf('Using complex differentials for gradient computation\n');
        end
        funObj = @autoGrad;
    else
        if debug
            fprintf('Using complex differentials for Hessian computation\n');
        end
        funObj = @autoHess;
    end
    
    if method == 5
        if debug
            fprintf('Turning off the use of complex differentials\n');
        end
        useComplex = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%% CONSTRAINT-RELEVANT CODE
% Form Projected Objective Function
varargin = {LB,UB,funObj,varargin{:}};
funObj = @projectedFunObj;

% Take projection of initial point
x = projectBounds(x0,LB,UB);

% Only backtracking line-searches can be used
LS = min(LS,2);
%%%%%%%%%%%%%%%%%%%%%%%

% Evaluate Initial Point
if method < 7
    [f,g] = funObj(x,varargin{:});
else
    [f,g,H] = funObj(x,varargin{:});
end
funEvals = 1;


% Output Log
if verboseI
    fprintf('%10s %10s %15s %15s %15s\n','Iteration','FunEvals','Step Length','Function Val','Opt Cond');
end

% Perform up to a maximum of 'maxIter' descent steps:
for i = 1:maxIter

    % Check Optimality Condition
    if sum(abs(g)) <= tolFun
        exitflag=1;
        msg = 'Optimality Condition below TolFun';
        break;
    end

%%%%%%%%%%%%%%%%%%%%%%% CONSTRAINT-RELEVANT CODE
    % Compute Free (working) Set
    free = ones(p,1);
    free((x < LB+1e-4) & g >= 0) = 0;
    free((x > UB-1e-4) & g <= 0) = 0;

    % Check that we have at least one free variable
    if sum(free==1) == 0
        exitflag=2;
        msg = 'No further progress possible, all variables active at upper/lower bound';
        break;
    end
%%%%%%%%%%%%%%%%%%%%%%%

    % ****************** COMPUTE DESCENT DIRECTION *****************

    d = zeros(p,1);
    if method == 0

        % Steepest Descent

        d(free==1) = -g(free==1);
    elseif method == 7

        % Newton search direction

        % Attempt to perform a Cholesky factorization of the Hessian
        [R,posDef] = chol(H(free==1,free==1));

        % If the Cholesky factorization was successful, then the Hessian is
        % positive definite, solve the system
        if posDef == 0
            d(free==1) = -R\(R'\g(free==1));

        else
            % otherwise, adjust the Hessian to be positive definite based on the
            % minimum eigenvalue, and solve with QR
            % (expensive, we don't want to do this very much)
            if debug
                fprintf('Adjusting Hessian\n');
            end
            H = H + eye(sum(free==1)) * max(0,1e-12 - min(real(eig(H(free==1,free==1)))));
            d(free==1) = -H\g(free==1);
        end
    elseif method == 3
        % L-BFGS

        % Update the direction and step sizes

        if i == 1
            d(free==1) = -g(free==1); % Initially use steepest descent direction
            old_dirs = zeros(p,0);
            old_stps = zeros(p,0);
            Hdiag = 1;
        else
            [old_dirs,old_stps,Hdiag] = lbfgsUpdate(g-g_old,x-x_old,corrections,debug,old_dirs,old_stps,Hdiag);
            curvSat = sum(old_dirs(free==1,:).*old_stps(free==1,:)) > 1e-10;
            d(free==1) = lbfgs(-g(free==1),old_dirs(free==1,curvSat),old_stps(free==1,curvSat),Hdiag);        
        end
        g_old = g;
        x_old = x;


    end

    % ****************** COMPUTE STEP LENGTH ************************

    % Directional Derivative
    f_old = f;
    gtd = g'*d;

    % Check that progress can be made along direction
    if gtd > -tolX
        exitflag=2;
        msg = 'Directional Derivative below TolX';
        break;
    end

    % Select Initial Guess

    if i == 1
        if method < 5
            t = min(1,1/sum(abs(g)));
        else
            t = 1;
        end
        f_old = f;
        gtd_old = gtd;
    else
        if LS_init == 0
            % Newton step
            t = 1;
        elseif LS_init == 1
            % Close to previous step length
            t = t*min(2,(gtd_old)/(gtd));
            gtd_old = gtd;
        else
            % Quadratic Initialization based on {f,g} and previous f
            t = min(1,2*(f-f_old)/(gtd));
            f_old = f;
        end

        if t <= 0
            t = 1;
        end
    end

    % Bound the initial step size
    if method < 5
        t = min(t,1e4/(1+sum(abs(g))));
    end

    % Line Search (Armijo Backtracking)

    % Perform Backtracking line search
    if method == 7
        [t,x,f,g,LSfunEvals,H] = ArmijoBacktrack(x,t,d,f,f,g,gtd,c1,LS,tolX,debug,doPlot,LS_saveHessianComp,funObj,varargin{:});
    else
        [t,x,f,g,LSfunEvals] = ArmijoBacktrack(x,t,d,f,f,g,gtd,c1,LS,tolX,debug,doPlot,1,funObj,varargin{:});
    end
    funEvals = funEvals + LSfunEvals;

    
    
%%%%%%%%%%%%%%%%%%%%%%% CONSTRAINT-RELEVANT CODE
    % Take projection of iterate
    x = projectBounds(x,LB,UB);
%%%%%%%%%%%%%%%%%%%%%%%

    % Output iteration information
    if verboseI
        fprintf('%10d %10d %15.5e %15.5e %15.5e\n',i,funEvals,t,f,sum(abs(g)));
    end

    % ******************* Check for lack of progress *******************

    if sum(abs(t*d)) <= tolX
        exitflag=2;
        msg = 'Step Size below TolX';
        break;
    end

    if abs(f-f_old) < tolFun
        exitflag=2;
        msg = 'Function Value changing by less than TolFun';
        break;
    end

    % ******** Check for going over iteration/evaluation limit *******************

    if funEvals > maxFunEvals
        exitflag = 0;
        msg = 'Exceeded Maximum Number of Function Evaluations';
        break;
    end

    if i == maxIter
        exitflag = 0;
        msg='Exceeded Maximum Number of Iterations';
    end

end

if verbose
    fprintf('%s\n',msg);
end
if nargout > 3
    output = struct('iterations',i,'funcCount',funEvals,...
        'algorithm',method,'firstorderopt',sum(abs(g)),'message',msg);
end

end

function [legal] = isLegal(v)
legal = sum(any(imag(v(:))))==0 & sum(isnan(v(:)))==0 & sum(isinf(v(:)))==0;
end

function [x] = projectBounds(x,LB,UB)
x(x < LB) = LB(x < LB);
x(x > UB) = UB(x > UB);
end

function [f,g,H] = projectedFunObj(x,LB,UB,funObj,varargin)

x = projectBounds(x,LB,UB);

if nargout == 1
    f = funObj(x,varargin{:});
elseif nargout == 2
    [f,g] = funObj(x,varargin{:});
else
    [f,g,H] = funObj(x,varargin{:});
end

end