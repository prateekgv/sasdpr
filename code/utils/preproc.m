function y = preproc(r, M, x)
% y = preproc(r, M, x)
%
% Preprocess signal x prior to SASS in order to reduce
% the potential occurrence of transients at start and end of signal.
%
% The preprocessing generates a new data vector of size N + 2M, where the
% 2M data points are obtained by the reflection over y-axis followed by
% reflection over x-axis
% 
% INPUT
%    r : degree of polynomial
%    M : number of end-point samples to use for fitting the polynomial
%    x : 1D signal 
%
% OUTPUT
%    y : 1D signal

% Ivan Selesnick,  2012


% convert to row vector

trans = false;
if size(x, 1) > 1
    x = x.';
    trans = true;
end

N = length(x);
% x1 = [nan(1,M), x, nan(1,M)];


p = polyfit(1:M, x(1:M), r);
y_start = polyval(p, -(M:-1:1));
p = polyfit(N-M+1:N, x(N-M+1:N), r);
y_end = polyval(p,N+1:N+M);

y = [y_start, x, y_end];
% plot(y); hold on; plot(x1);

if trans
    y = y.';
end

