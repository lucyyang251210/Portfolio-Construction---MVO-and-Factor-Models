function  [mu, Q] = BSS(returns, factRet, lambda, K)
    
    % Use this function for the BSS model. Note that you will not use 
    % lambda in this model (lambda is for LASSO).
    %
    % You should use an optimizer to solve this problem. Be sure to comment 
    % on your code to (briefly) explain your procedure.
    
 
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    factRet = factRet;
    returns = returns;
    
    % Number of factors
    n = size(factRet,2);
    

    % Number of observations;
    N = size(factRet, 1);
    % Calculate the factor expected excess return from historical data using
    % the geometric mean
    %mu = (geomean(factRet + 1) - 1)';
    % Calculate the asset covariance matrix
    %Q = cov(factRet);
    %targetRet = geomean(factRet + 1) - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 2: Cardinality and buy-in thresholds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total of k out of 8 total factor
    k = K;    
%-------------------------------------------------------------------------- 
% 3.1 Inequality constraints:
% Gurobi accepts inequality constraints of the form "A x <= b" and 
% "A x >= b". However, for consistency, we will keep all constraints as 
% "A x <= b"
%--------------------------------------------------------------------------

    alphaPlace = ones(N,1);
    X = [alphaPlace factRet];
    Q= X'* X;
    I = eye(n+1);    
    A = [-I 0*I;
        -I I];
    Q = [Q zeros(n+1);
        zeros(n+1,2*(n+1))];

% We must also define the right-hand side coefficients of the inequality 
% constraints, b:
    
    b = zeros(2*n+2,1);
    %display(targetRet);

%--------------------------------------------------------------------------
% 3.2 Equality constraints: 
% We will define our cardinality constraint as an equality (although this
% is not always the case, we could also define it as an inequality)
%--------------------------------------------------------------------------

% We only have 2 equality constraints: the cardinality constraint (sum of
% y's) and the budget constraint (sum of x's):
    Aeq = [zeros(1,n+1) ones(1,n+1)];

% We must also define the right-hand side coefficients of the equality 
% constraints, beq:
    beq = [k];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PART 4: Setup Gurobi model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% 4.1 Define the model variables and assign them a name
%--------------------------------------------------------------------------

% Define the variable types:'C' defines a continuous variable, 'B' defines
% a binary variable
    varTypes = [repmat('C', n+1, 1); repmat('C', n+1, 1)];

% Input the lower and upper bounds. Since our lower buy-in threshold is
% 0.05, this means we are not allowed to short-sell.
    lb = zeros(2*n+2, 1);
    ub = ones(2*n+2, 1); 

% Assign tags for the model variables.

% Append '_c' to cont var. names
    %namesCont = cellfun(@(c)[c '_c'], tickers(1:n), 'uni', false); 

% Append '_b' to binary var. names
    %namesBin = cellfun(@(c)[c '_b'], tickers(1:n), 'uni', false); 

% Combine both name vectors
    %names = [namesCont namesBin];

%--------------------------------------------------------------------------
% 4.1 Setup the Gurobi model
%--------------------------------------------------------------------------
    clear model;

% Assign the variable names
    %model.varnames = names;

% Gurobi accepts an objective function of the following form:
% f(x) = x' Q x + c' x 

% Define the Q matrix in the objective 
    model.Q = sparse(Q);

% define the c vector in the objective (which is a vector of zeros since
% there is no linear term in our objective)
    model.obj = zeros(1, 2*n+2);

% Gurobi only accepts a single A matrix, with both inequality and equality
% constraints
    model.A = [sparse(A); sparse(Aeq)];

% Define the right-hand side vector b
    model.rhs = full([b; beq]);

% Indicate whether the constraints are ">=", "<=", or "="
    model.sense = [ repmat('<', (2*n+2), 1) ; repmat('=', 1, 1) ];

% Define the variable type (continuous, integer, or binary)
    model.vtype = varTypes;

% Define the variable upper and lower bounds
    model.lb = lb;
    model.ub = ub;

% Set some Gurobi parameters to limit the runtime and to avoid printing the
% output to the console. 
    clear params;
    params.TimeLimit = 100;
    params.OutputFlag = 0;
 
    results = gurobi(model,params);
    
    mu = results.objval + [(-2) * returns' * X]';
    Q = sparse(Q);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End


