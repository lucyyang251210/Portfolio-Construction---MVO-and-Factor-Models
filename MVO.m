function  x = MVO(mu, Q, targetRet)
    
    % Use this function to construct your MVO portfolio subject to the
    % target return, with short sales disallowed. 
    %
    % You may use quadprog, Gurobi, or any other optimizer you are familiar
    % with. Just be sure to comment on your code to (briefly) explain your
    % procedure.

    % Find the total number of assets
    n = size(Q,1); 

    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    H = Q;
    f = [];
    A = - mu.';
    b = [-targetRet];
    Aeq = ones(1,n);
    beq = [1];
    lb = zeros(1,n);
    ub = ones(1,n);
    x =  quadprog(H, f, A, b, Aeq, beq, lb, ub); % Optimal asset weights
    %----------------------------------------------------------------------
    
end