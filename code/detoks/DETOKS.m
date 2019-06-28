function [x,s,f,cost] = DETOKS(y,fs,d,fc,lam1,lam2,lam3,Nit,mu)
% function [x,s,f,cost] = DETOKS(y,fs,d,fc,lam1,lam2,lam3,Nit,mu)
%Spindle detection using sparse regularization
%
% input - 
%         y - 'C3-A1' channel
%         fs - sampling frequency 
%         d - degree of the low pass filter
%         fc - cut-off frequency of the low pass filter
%         lam1,lam2,lam3 - regularization parameters
%         mu - ADMM convergence parameter
%         Nit - Number of iterations
%         `windowsize for stft is 1 sec'
%         
% output - 
%         x - sparse signal with sparse first order derivative
%         s - spindle component
%         f - low frequency component
%         
%   Please cite as: 
%   A. Parekh, I.W. Selesnick, D. Rapaport, I.Ayappa, Detection of
%   K-complexes and sleep spindles (DETOKS) using sparse optimization,
%   Journal of Neuroscience Methods, vol. 251, pp. 37-46, 2015. 
%   
%   Ankit Parekh(ankit.parekh@nyu.edu)
%   Copyright (c) 2015. 


y = y(:);                               % Convert to column vector
cost = zeros(Nit,1);                    % Cost function history
N = length(y);   
[A, B] = ABfilt(d, fc, N);              % Create the low-pass filter
H = @(x) A\(B*x);                       % H: high-pass filter
G = mu*(A*A') + 2*(B*B');               % G: Banded system

% For STFT parameters see I. W. Selesnick, “Short-time Fourier transform and its inverse,” 
% http://cnx.org/content/m32294/, 2009.

    [A1,A1H,~] = MakeTransforms('STFT',N, [2^nextpow2(fs) 4 2 2^nextpow2(fs)]); 

c = A1H(y);
x = zeros(N,1);
d1 = zeros(N,1);
d2 = A1H(y);

b = (1/mu) * B' * ((A*A')\(B*y));
Ab = A1H(b);

for i = 1:Nit
    g1 = b + (x+d1);
    g2 = Ab + (c+d2);
    u1 = g1 - B' * (G\(B*(g1 + A1(g2)')));
    u2 = g2 - A1H(B' * (G\(B*(g1 + A1(g2)'))));
    x = soft(tvd(u1-d1,N,lam2/(mu)),lam1/(mu))';
    c = soft(u2-d2,lam3/(mu));
    d1 = d1 - (u1 - x);
    d2 = d2 - (u2-c);
    cost(i) = 0.5*sum(abs(H(y-x-real(A1(c)'))).^2) +...
              lam1 * sum(abs(x)) + lam2 * sum(abs(diff(x))) + ... 
              lam3 * sum(abs(c(:))) ;
end

s = real(A1(c)');
f = y - x - s - H(y-x-s);

% Bandpsas filtering
[B0,A0] = butter(4, [11.5/(fs/2) 15.5/(fs/2)]);
s = filtfilt(B0,A0,s);

% Bandpass filtering
% [B0,A0] = butter(4, [0.6/(fs/2) 1.4/(fs/2)]);
% f = filtfilt(B0,A0,f);