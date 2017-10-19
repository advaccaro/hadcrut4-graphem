function [bias_sq, variance, MSE] = MSE_decomp(X_target, X)
% function [bias_sq, variance, MSE] = MSE_decomp(X_target, X)
%      Bias-Variance Decomposition of a matrix X
%    Written by Dominique Guillot (USC), Feb 2011
%    Edited by Julien EMile-Geay (USC), June 2011
% =========================================================

MSE = mean(mean((X - X_target).^2));

bias = mean(mean(X-X_target));

bias_sq = bias^2;

variance = MSE - bias_sq;  % Using Var(X) = E(X^2) - E^2(X).

end
