function  [mu, Q] = LASSO(returns, factRet, lambda, K)
    
    % Use this function for the LASSO model. Note that you will not use K 
    % in this model (K is for BSS).
    %
    % You should use an optimizer to solve this problem. Be sure to comment 
    % on your code to (briefly) explain your procedure.
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    factRet = factRet;
    returns = returns;
    % rowNum is the number of observation
    [rowNum, assetNum] = size(returns);
    
    % m is rowNum; 8 is facNum
    [rowNum, facNum] = size(factRet);
    
    % we can write the problem to below format
    % min B'X'XB + c'B = B'X'XB + (- 2) r'X B + lambdaVec* B
    % s.t. B unrestircitved A should be indentity matrix here
    % let B = B+ - B- B+ B- >=0
    
    t = rowNum;
    n = assetNum;
    p = facNum;
    alphaPlace = ones(rowNum,1);
    X = [alphaPlace factRet]; % m * p+1 matrix
    r = returns; % m * n matrix
    
    % as Bi is set to B+ -B- B should be 2p + 2 
    
    H = X'* X;
    % H need to be a 2p+2 matrix ??? 
    H = [H zeros(p+1);
        zeros(p+1,2*(p+1))];
    
    % Now we want to get A : xi+ + xi- < = s xi+ >=0
    I = eye(p+1);
    
    A = [I -I;
        -I -I];
    b = zeros(2*p +2,1);
    
    lambdaVec = ones(1, 1+p) * lambda;
    
    % use quodraprog to solve for each asset i : n
    % min 1/2 B'HB + c'B st AB <= 0
    % equavalent B'Hb + 2c'B st Ab<= 0 b+,b- >=0
    %betaArray = zeros(n, p+1) ;
    for ii = 1:n
        ri = r(:,ii);
        f = [(-2) * ri' * X  lambdaVec]';
        beta= quadprog(H,f,A,b,[],[],[],[]);
        betaArray(:,ii) = beta;
    end
    
    alpha = betaArray(1,:);
    B = betaArray(2:facNum+1,:);
    meanFac = mean(factRet); % 1 by 8 for 8 factors 

    E = factRet * B - returns;
    covE = 1/(rowNum-facNum-1)* E' * E;
    % mu =          % n x 1 vector of asset exp. returns
    mu = alpha +  meanFac * B;
    mu = mu';
    % Q  =          % n x n asset covariance matrix
    covFac = cov(factRet);
    D = diag(covE);
    Q = B'*covFac*B + D;
    %----------------------------------------------------------------------
    
end