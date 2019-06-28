function [x, u, cost, v] = sass_L1(y, d, fc, K, lam, u_init)
% x = sass_L1(y, d, fc, K, lam)
%
% Sparsity-assisted signal smoothing (SASS)
% using L1 norm as sparsity penalty.
%
% Signal model: y = f + g + noise
%   f : low-pass signal
%   g : sparse order-K derivative
%
% INPUT
%   y - noisy data
%   d - filter degree parameter (d = 1, 2, 3)
%   fc - cut-off frequency (normalized, 0 < fc < 0.5)
%   K - order of sparse derivative (1 <= K <= 2d)
%   lam - regularization parameter
%
% OUTPUT
%   x - output of SASS algorithm
%
% Use [x, cost, u, v] = sass_L1(...) to obtain
%   u - sparse signal
%   cost - cost function history
%   v - filtered sparse signal, A\(B1*u)

% Ivan Selesnick,  NYU Tandon School of Engineering
% 2012, revised 2016
% References:
% [1] Sparsity-Assisted Signal Smoothing
% [2] Sparsity-Assisted Signal Smoothing (Revisited)
% http://eeweb.poly.edu/iselesni/sass

MAX_ITER = 10000;
TOL_STOP = 1e-4;

y = y(:);                           % Convert to column
cost = zeros(1, MAX_ITER);          % Cost function history
N = length(y);

[A, B, B1, D] = ABfilt(d, fc, N, K); % Banded filter matrices

H = @(x) A\(B*x);                   % H : high-pass filter
AAT = A*A';                         % A*A' : banded matrix [sparse]
Hy = H(y);
b = B1'*(AAT\(B*y));
if exist('u_init', 'var')           % initialization
    u = u_init;
else
    % u = diff(y, K);
    u = D * y;
end

L = size(B1, 2);

iter = 0;
old_u = u;
last_iter = false;

while not(last_iter)
    
    iter = iter + 1;
    
    Lam = spdiags(abs(u)/lam, 0, L, L);       % Lam : diagonal matrix
    Q = AAT + B1*Lam*B1';                     % Q : banded matrix
    u = Lam * (b - (B1'*(Q\(B1*(Lam*b)))));   % Update
    cost(iter) = 0.5 * sum( abs(Hy-A\(B1*u)).^2 ) + lam * sum(abs(u));

    delta_u = max(abs( u - old_u )) / (max(abs(old_u))+eps);
    old_u = u;
    
    if (delta_u <= TOL_STOP) || (iter >= MAX_ITER)
        last_iter = true;
    end
    
end

v = A\(B1*u);
x = y - Hy + v;           % Total solution L(y) + A\(B1*u)

% Check for zero-locking
r = B1'*(AAT\(B*y-B1*u)) / lam;
k = abs(r) > 1.02;               % Use margin of 0.02
NZL = sum(k);                    % NZL : Number of falsely locked zeros
% fprintf('%d incorrectly locked zeros detected.\n', NZL);

cost = cost(1:iter);
