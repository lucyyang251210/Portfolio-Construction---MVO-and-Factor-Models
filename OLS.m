function  [mu, Q] = OLS(returns, factRet, lambda, K)
    
    % Use this function to perform an OLS regression. Note that you will 
    % not use lambda or K in thismodel (lambda is for LASSO, and K is for
    % BSS).
 
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % y = X beta + epsilon
    % X is the factors return

    factReturn  = factRet;
    % y is the returns that need to estimate it's a m by 20 matrix
    %y = returns{:,:};
    y = returns;
    % m is rowNum i.e. # of observation; 20 is assetNum
    [rowNum, assetNum] = size(returns);
    
    % m is rowNum; 8 is facNum
    [rowNum, facNum] = size(factReturn);
    % y = alpha + X beta + epsilon need to incorporate alpha
    
    alphaPlace = ones(rowNum,1);
    % add a column to X so X is a m by 8+1=9 matrix 
    
    X = [alphaPlace factReturn];
    betaArray = (X' * X) \ X' * y;
    
    alpha = betaArray(1,:);
    B = betaArray(2:facNum+1,:);
    
    % m by 9 * 9 by 20 = m by 20
    E = X * betaArray - y;
    % Take a square of each element in residuals
    
    %E_sqr=E.^2;
    % Compute residuel variance
    %Cov_E = 1/(Num_Observation-Num_factor-1)*sum(E_sqr(1:end,:));
    sqrE = E.^2 ;
    % Compute residuel variance
    covE = 1/(rowNum-facNum-1)* sum(sqrE(1:end,:));

      
    % mu =          % n x 1 vector of asset exp. returns
    meanFac = mean(factReturn); % 1 by 8 for 8 factors 
    mu = alpha +  meanFac * B;
    mu = mu';
    
    % Q  =          % n x n asset covariance matrix
    % asset cov = beta * factor cov * beta' + D
    covFac = cov(factReturn);
    % D is diagonal . it is the residual variances's diagonal elements
    D = diag(covE);
    Q = B' * covFac * B + D;
    
    %B
    %----------------------------------------------------------------------
    
end