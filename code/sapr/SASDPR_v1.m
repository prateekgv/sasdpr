%% function [x,s,f,c,cost] = SASDPR_v1(y,fs,fc,lam1,lam2,lam3,Nit,mu)
%
% Inputs:
%   y       - input signal
%   fs      - sampling rate
%   fc      - cutoff freq. of high-pass filter
%   lam1    - regularization parameter of the sparse vector
%   lam2    - regularization parameter of the sparse derivative vector
%   lam3    - regularization parameter of the oscillatory signal
%   Nit     - maximum number of iterations
%   mu      - rate of convergence parameter
%
% Outputs:
%   x       - sparse and sparse derivative signal
%   s       - oscillatory pattern
%   f       - low-frequency signal
%   c       - short time Fourier transform dictionary
%   cost    - cost function
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%


function [x,s,f,c,cost] = SASDPR_v1(y,fs,fc,lam1,lam2,lam3,Nit,mu)

y = y(:);                                                                   % Convert to column vector
N = length(y);
cost = zeros(Nit,1);

NB = 4; NA = 4; type = 'high'; K = 0;
wn = 0.1; wc = fc; deg = [NB,NA];                                                             
[H_hp,~,~,~,~,~,~] = IIR_ABfilt(deg, N, [wn,wc], type);

% band-pass filter
NB = 4; NA = 4; type = 'band'; K = 0;
wn = 4/(fs/2); wc = 13/(fs/2); deg = [NB,NA];
[B_bp,~,~,~,~,~,~] = IIR_ABfilt(deg, N, [wn,wc], type);

H_zp = H_hp'*H_hp;                                                          % H_zp: zero-phase high-pass filter
B_zp = B_bp'*B_bp;                                                          % B_zp: zero-phase band-pass filter
Hp = @(x) (H_zp*x);                                                         % Hp: high-pass filter operator

[A1,A1H,~] = MakeTransforms('STFT',N, [2^nextpow2(fs) 4 2 2^nextpow2(fs)]);


x = H_zp*y;
c = A1H(B_zp*y);

d1 = zeros(N,1);
d2 = zeros(size(c));

G = H_zp'*H_zp + B_zp'*B_zp + (mu)*eye(size(B_zp,1));
F = inv(G);
BF = B_zp*F;
HF = H_zp*F;

b1 = (1/mu)*(H_zp'*H_zp)*y;
b2 = A1H(B_zp*b1);

for i = 1: Nit
    
    g1 = x + d1 + b1;
    g2 = c + d2 + b2;
    
    g = B_zp*(A1(g2)') + H_zp*g1;
    u1 = g1 - HF*g;
    u2 = g2 - A1H(BF*g);
    
    x = soft(tvd(u1-d1,N,lam2/(mu)),lam1/(mu))';
    c = soft(u2-d2,lam3/(mu));
    
    d1 = d1 - (u1 - x);
    d2 = d2 - (u2 - c);
    
    cost(i) = 0.5*sum(abs(H_zp*(y-x-B_zp*real(A1(c)'))).^2) +...
              lam1 * sum(abs(x)) + lam2 * sum(abs(diff(x))) + ... 
              lam3 * sum(abs(c(:))) ;
end

s = B_zp*real(A1(c)');
f = y - x - s - Hp(y-x-s);

end

