%%Gaussian Maximum Likelihood Estimator

function [Mu, Sigma] = GaussianMLEstimator(X)
A = sum(X)/length(X);                   %for multivariate case, the estimate value of mu is
b = 0;                                  %the mean vector of Xk
for i = 1:length(X)
    a = (X(i,:) - A)'*(X(i,:) - A);
    b = b + a;                          %the sum of each covariance matrix
end
B = b/length(X);                        %devided by n then get the square root
Mu = A;
Sigma = B;
