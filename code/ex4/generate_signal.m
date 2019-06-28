function [y, f, s, x, w] = generate_signal( fs, sigma )

% Initialize noise level and signal length
rng('default')
N = 10*fs;
n = 0:N-1;

% Generate synthetic components
% low-frequency
f = 0.1;
f = sin(2*f/fs*pi*n);      

% oscillatory
s = zeros(size(n));
o = 13;
s(200+(1:fs)) = sin(2*pi*o/fs*(1:fs)) .* hamming(fs)';

% sparse piece-wise constant
x = zeros(size(n));
x(100:110) = -1;
x(400:420) = 1;

% add noise
w = sigma*randn(size(n));
y = f+s+x+w;

end

