%% function [kcomp,r_kcomp,cost] = SAPR_v1(y,fs,fc,lam0,lam1,mu,eta,Nit)
%
% Inputs:
%   y       - input signal
%   fs      - sampling rate
%   fc      - cutoff freq. of high-pass filter
%   lam0    - regularization parameter sparsity of wavelet coefficients
%   lam1    - regularization parameter first order difference
%   mu      - rate of convergence parameter
%   eta     - rate of convergence parameter
%   Nit     - maximum number of iterations
%
% Outputs:
%   kcomp   - band pass filtered detected patterns
%   r_kcomp - reconstructed patterns without band pass filtering
%   cost    - cost function
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%
function [kcomp,r_kcomp,cost] = SAPR_v1(y,fs,fc,lam0,lam1,mu,eta,Nit)

y = y(:);
cost = zeros(Nit,1);                    % Cost function history
N = length(y);

% High-Pass Filter at 0.6 Hz
NB = 4; NA = 4; type = 'high'; K = 0;
wn = 0.1; wc = fc; deg = [NB,NA];                                                             
[H_hp,~,~,~,~,~,~] = IIR_ABfilt(deg, N, [wn,wc], type);

% Band-Pass Filter (0.6-1.8 Hz) for K-Complex
NB = 4; NA = 4; type = 'band'; K = 0;
wn = 1.4/(fs/2); wc = 1.3/(fs/2); deg = [NB,NA];
[B_bp,~,~,~,~,~,~] = IIR_ABfilt(deg, N, [wn,wc], type);

H_zp = H_hp'*H_hp;                                                          % H_zp: zero-phase high-pass filter
B_zp = B_bp'*B_bp;                                                          % B_zp: zero-phase band-pass filter
Hp = @(x) (H_zp*x);                                                         % Hp: high-pass filter operator

G = B_zp'*B_zp + (mu)*eye(size(B_zp,1));
BF = B_zp/G;

[A1,A1H,~] = MakeTransforms('DWT2',N, [2^nextpow2(fs) 4 4 2^nextpow2(fs)]); 

k = A1H(B_zp*y);
v = k;

d1 = zeros(size(k));
d3 = zeros(size(v));

b = (1/mu)*H_zp*y;
b1 = A1H(B_zp*b);


for i = 1:Nit
    
    g1 = b1 + d1 + k;                           % update g1
    g = B_zp*(A1(g1)');

    temp1 = BF*g;
    u1 = g1 - A1H(temp1);                       % update u1
    
    p = (mu*(u1-d1)+eta*(v-d3))/(mu+eta);       % update p
    k = soft(p,lam0/(mu+eta));                  % update k
    
    m = d3 + k;                                 % update m
    temp3 = A1(m)';
    temp4 = tvd(temp3,N,lam1/eta);
    temp5 = temp4' - A1(m)';
    v = m + A1H(temp5);                         % update v
    
    d1 = d1 - (u1 - k);                         % update d1
    d3 = d3 - (v - k);                          % update d3
    
    cost(i) = 0.5*sum(abs(Hp(y) - B_zp*(A1(k))' )).^2 +...
              lam0 * sum(abs(k(:))) + lam1 * sum(abs(diff(A1(k))));
end

r_kcomp = real(A1(k)');
kcomp = B_zp*real(A1(k)');

end