%% function [z,cost] = FISTA_L1( A, b, lambda, L, positive, tolerance, max_iter, verbose)
% Description:
%             FISTA_LASSO_solve Solves LASSO via FISTA.
%             Solves min_x 1/2||Ax-b||_2^2 + \lambda ||x||_1
%
% Inputs:
%   A       - Input matrix
%   b       - Input vector
%   lambda  - Regularization parameter
%   L       - Lipschitz constant
%
% Outputs:
%   z       - Output sparse vector
%   cost    - Output cost function
%
% Author: Mianzhi Wang
% Modified by: Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%
function [z,cost] = FISTA_L1( A, b, lambda, L, positive, tolerance, max_iter, verbose)
%FISTA_LASSO_solve Solves LASSO via FISTA.
% Solves min_x 1/2||Ax-b||_2^2 + \lambda ||x||_1

if nargin <= 7
    verbose = false;
end
if nargin == 3
    L = [];
    positive = false;
    tolerance = 1e-4;
    max_iter = 50;
elseif nargin == 4
    positive = false;
    tolerance = 1e-4;
    max_iter = 50;
elseif nargin == 5
    tolerance = 1e-4;
    max_iter = 50;
elseif nargin == 6
    max_iter = 50;
end

if ~isreal(A) || ~isreal(b)
    error('A and y must be real.');
end

if isempty(L)
    L = eigs(A*A', 1);
end

cost = zeros(1,max_iter);

t = 1;
n = size(A,2);
% x = pinv(A) * b;
x = zeros(size(A,2), 1);
y = x;
k = 0;
ATA = A' * A;
r = A' * b;
mu = lambda / L;
eta = 1.2;
last_obj_func_val = 0.5*sumsqr(A*x - b) + lambda*norm(x, 1);
while true

    x_old = x;
    % proximal + shrink & threshold
    p = y - (ATA * y - r) / L;
    if positive
        x = soft_thresholding_p(p, mu);
    else
    	x = soft_thresholding(p, mu);
    end
    if k > 0 && norm(x - x_old) / norm(x_old) < tolerance
        break;
    end
    
    t_old = t;
    t = (1 + sqrt(1 + 4*t*t)) / 2;
    y = x + (t_old - 1) * (x - x_old) / t;
    
    k = k + 1;
    fit_val = 0.5*sumsqr(A*x - b);
    obj_func_val = fit_val + lambda*norm(x, 1);
    cost(1,k) = obj_func_val;
    improvement = abs(obj_func_val - last_obj_func_val) / last_obj_func_val;
    if verbose
        fprintf('%s obj = %.6e, R = %.6e, L = %.6e, I = %.6e \n', ...
                sprintf('[%d]', k), ...
                obj_func_val, fit_val, L, improvement);
    end
    last_obj_func_val = obj_func_val;
    
    if k > max_iter || improvement < tolerance
        break;
    end
    
end

z = x;
cost = cost(2:k);

end


function s = soft_thresholding(x, t)
n = length(x);
s = zeros(n, 1);
for i = 1:n
    if x(i) > t
        s(i) = x(i) - t;
    elseif x(i) < -t
        s(i) = x(i) + t;
    end
end
end

function s = soft_thresholding_p(x, t)
n = length(x);
s = zeros(n, 1);
for i = 1:n
    if x(i) > t
        s(i) = x(i) - t;
    end
end
end