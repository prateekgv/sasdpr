%% function [x, x1, x2, cost, g] = sasd_L1(y, K, lam, H, H1)
%
% Inputs:
%   y       - noisy signal
%   K       - sparsity or derivative order
%   lam     - regularization parameter
%   H       - zero-phase high pass filter
%   H1      - K-order sparse derivative of the high-pass filter
%
% Outputs:
%   x1      - reconstructed low-frequency signal
%   x2      - reconstructed k-order sparse derivative
%   x       - reconstructed signal
%   cost    - cost function
%   g       - optimization plot
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%


function [x, x1, x2, cost, g] = sasd_L1(y, K, lam, H, H1)

% x - denoised signal
% v - sparse signal
% x1 - reconstructed low-frequency signal
% x2 - reconstructed k-order sparse derivative

y = y(:);  
N = length(y);

L = @(x) x - H'*H*x;

HTH1 = full(H'*H1);
Hy = full(H'*H*y);
[v,cost] = FISTA_L1(HTH1, Hy, lam, [], false, 1e-6, 10000);

S = ones(N,N-K);
S = tril(S,-K);
x2 = S*v;
x1 = L(y - x2);

v1 = H'*(H1*v);
x = L(y) + v1;           % Total solution L(y) + v1

g = (H1'*H)*(Hy - H'*H1*v);

end