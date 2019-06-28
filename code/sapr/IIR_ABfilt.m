%==========================================================================
% function [H,b_tf,a_tf] = IIR_ABfilt(NB,NA,N,omega,type)
%==========================================================================
% @author     : Prateek Gundannavar
% @email      : prateek@wustl.edu
% @description: This function develops filters as matrices. First, we
%               develop a prototype low-pass digital genearlized
%               Butterworth filter and convert it into a composite filter
%               depeneding on the type of the filter. Possible types
%               include low-pass, high-pass, and band-pass filters
% @inputs     :
%             - deg     [NA,NB] contains the deg of num and denom
%                       polynomial
%             - N       The size of the matrix filter
%             - omega   omega(1) - The cut-off frequency of the prototype
%                       LPF or equivalently the band-width of the BPF
%                       omega(2) - The cut-off frequency of the composite
%                       LPF/HPF or the central frequency of the band-pass
%             - type    'low'  - Composite LPF
%                       'high' - Composite HPF
%                       'band' - Composite BPF
%             - K       Degree of the sparse derivative value
% @outputs    :
%             - H       Lower-triangle Hankel matrix consisting of impulse
%                       response of the composite filter G(z)
%             - H1      Factorization matrix of size N x N-K
%             - D       K-sparse derivative matrix of size N-K x N
%             - b_tf    The numerator coefficients of G(z)
%             - a_tf    The denominator coefficients of G(z)
%             - err     The error obtained after factorization
%==========================================================================
function [H,H1,D,b_tf,a_tf,H2norm,err] = IIR_ABfilt(deg,N,omega,type,K,verbose)

%% Check the sparse derivative value
if nargin < 5
    K = 0;
    err = 0;
    H2norm = [];
end

if nargin < 6
    verbose = false;
end

%% Generate a prototype low-pass maxflat filter
if verbose
    disp('Creating a prototype maxflat low-pass filter H(z)');
end
[b_main,a_main,b1,b2] = maxflat(deg(2),deg(1),omega(1));

% The sparse derivative degree cannot exceed the number of zeros at z=-1.
if (K > length(b1)-1)
    error('The max. sparse derivative degree possible is %d',length(b1)-1);
end

%% State-space represntation of the prototype filter
if (deg(2) <= deg(1))
    [sos,g] = tf2sos(b_main,a_main);
    [A_main,B_main,C_main,D_main] = sos2ss(sos,g);
else
    a_main = [a_main,zeros(1,deg(2)-deg(1))];
    [A_main,B_main,C_main,D_main] = tf2ss(b_main,a_main);
end

%% Apply balance transformation
if verbose
    disp('Appling balance transformation to the prototype filter H(z)');
end
[A_bal,B_bal,C_bal,D_bal,~] = balancing_transformation(A_main,B_main,C_main,D_main);

%% Generate the composite filter
[A_tf,B_tf,C_tf,D_tf] = generate_composite_filter(A_bal,B_bal,C_bal,D_bal,omega,type);
[b_tf,a_tf] = ss2tf(A_tf,B_tf,C_tf,D_tf);

%% Generate the impulse response matrix
if verbose
    disp('Generating an impulse response matrix for G(z)');
end
H = generate_filt_matrix_zero_inital(A_tf,B_tf,C_tf,D_tf,N,N);

%% Create Factorization Matrix D (Only for high-pass and band-pass filters)
[D,d] = generate_k_sparse_derivative_matrix(K,N);                           % N-K x N (derivative matrix)
if ~strcmp(type,'low') && K~=0

    [b_aux,r_aux] = deconv(b_tf,d);
    if r_aux(end) < 1e-6
        [A_aux,B_aux,C_aux,D_aux] = tf2ss(b_aux,a_tf);
        [A_bal,B_bal,C_bal,D_bal] = balancing_transformation(A_aux,B_aux,C_aux,D_aux);
        H1_init = generate_filt_matrix_zero_inital(A_bal,B_bal,C_bal,D_bal,N,N-K);
        [H1,err] = acc_projected_grad_descent(H,D,N,K,H1_init);
        
        imp = zeros(N-K,1);
        imp(round(N/2)) = 1;                                                % imp : impulse signal (located at center to avoid transients)
        h2 = H'*H1*imp;
        H2norm = sqrt( sum( abs( h2 ).^2 ) );                               % H2norm
    else
        error('Factorization cannot be done as the residue is large');
    end
    
else
    H1 = [];
end

end


%% AUXILLARY FUNCTIONS

%==========================================================================
% function [H1,err] = acc_projected_grad_descent(H,D,N,K,H0)
%==========================================================================
% @description: Solves a quadratic program to obtain a lower-triangular
%               matrix which gives the matrix factorization of the filter
%
% @inputs     :
%             - H       N x N matrix representing the composite filter
%             - D       N-K x N sparse derivative matrix
%             - N       N size of the input vector
%             - K       K sparse-derivative value
%             - H0      Initial value of H1
% @outputs    :
%             - H1      N x N-K factorization matrix
%==========================================================================
function [H1,err] = acc_projected_grad_descent(H,D,N,K,H0,verbose)


if nargin <= 5
    verbose = false;
end

tolerance = 1e-4;
max_iter = 10000;

L = [];

HHT = lt_inner_prod(H');
HTH = lt_inner_prod(H);
DDT = D*D';
H1 = H0;

if isempty(L)
    [~,D1,~] = svd(HHT,'econ');
    [~,D2,~] = svd(full(DDT),'econ');
    L = D1(1,1)*D2(1,1);
end

t = 1;
x = L*vec(H1);
y = x;
k = 0;
r = vec(HHT*H*D');

last_obj_func_val = 0.5*norm(HTH - H'*H1*D,'fro');

while true
    
    x_old = x;
    p = (y - (vec(HHT*H1*DDT) - r )/L);      % proximal gradient
    H1 = reshape(p,N,N-K);                  
    H1 = tril(H1);                          % projection
    x = vec(H1);
    
    if k > 0 && norm(x - x_old) / norm(x_old) < tolerance
        break;
    end
    
    t_old = t;
    t = (1 + sqrt(1 + 4*t*t)) / 2;
    y = x + (t_old - 1) * (x - x_old) / t;
    
    k = k + 1;
    fit_val = 0.5*norm(HTH - H'*H1*D,'fro');
    obj_func_val = fit_val;
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

err = fit_val;


end

%==========================================================================
% D = generate_k_sparse_derivative_matrix(K,N)
%==========================================================================
% @description: Creates a sparse-derivative matrix D of order K
%
% @inputs     -
%             K - degree of sparse derivative
%             N - size of the input
% @outputs    -
%             D - N-K x N sparse matrix
%==========================================================================
function [D,d] = generate_k_sparse_derivative_matrix(K,N)
d = 1;
for i = 1:K
    d = conv(d, [-1 1]);
end
D = spdiags(d(ones(N,1), :), 0:K, N-K, N);               % D: banded matrix
end

%==========================================================================
% function [A_b,B_b,C_b,D_b,T] = balancing_transformation(A,B,C,D)
%==========================================================================
% @description: Applies balancing transformation to a given state-space
%               representation of a digital filter
%
% @inputs     -
%               A, B, C, D - State-space representation of the filter 
%               using tf2ss
%
% @outputs    -
%               A_b, B_b, C_b, D_b : Gramian preserving transformations
%               T                  : Transformation matrix
%==========================================================================
function [A_b,B_b,C_b,D_b,T] = balancing_transformation( A,B,C,D )

% Generate Gramian matrices
m = size(A,1);
K = zeros(m,m); W = zeros(m,m);
A_old = eye(size(A,1));
for i = 1:2000
    K = K + (A_old*B)*(A_old*B)';                   % Reachability
    W = W + (C*A_old)'*(C*A_old);               % Observability
    A_old = A*A_old;
end

[L_c,p1] = chol(K,'lower');
[L_o,p2] = chol(W,'lower');

% Perform Gramian transformation only when positive definite
if (p1 == 0) && (p2 == 0)
    [~,Sigma,V] = svd(L_o'*L_c);
    
    % Obtain the similarity transform
    T = L_c*V*inv(sqrt(Sigma));
    
    A_b = inv(T)*A*T;
    B_b = inv(T)*B;
    C_b = C*T;
    D_b = D;

else
    error('Gramian matrices are unstable.')
end


end


%==========================================================================
% function [A_tf,B_tf,C_tf,D_tf] = generate_composite_filter(A_bal,B_bal,...
%    C_bal,D_bal,omega,type)
%==========================================================================
%
% @description: Generate the composite filter based on the prototype
%               low-pass filter
%
% @inputs     :
%               - A_bal, B_bal, C_bal, D_bal - State-space representation 
%               of the prototype balanced low-pass filter
%               - omega - contains the cut-off frequency of the high-pass
%               or low-pass filter, and the central frequency of the
%               band-pass filter.
%               - type - type of the composite filter G(z)
%
% @outputs    :
%               A_tf, B_tf, C_tf, D_tf - Gramian preserving transformations
%               of the composite filter G(z)
%==========================================================================
function [A_tf,B_tf,C_tf,D_tf] = generate_composite_filter(A_bal,B_bal,...
    C_bal,D_bal,omega,type,verbose)

if nargin <= 6
    verbose = false;
end

% Determine the type of spectral transformation
switch type
    case 'low'
        wo = omega(1); wt = omega(2);
        eta = + sin(pi*(wo-wt)/2) / sin(pi*(wo+wt)/2);
        b_ap = [-eta 1];
        a_ap = [1 -eta];
        [A,B,C,D] = tf2ss(b_ap,a_ap);
        if verbose
            disp('Appling balance transformation to the all-pass filter 1/F(z)');
        end
        [alpha,beta,gamma,delta] = balancing_transformation(A,B,C,D);
    case 'high'
        wo = omega(1); wt = omega(2);
        eta = - cos(pi*(wo+wt)/2) / cos(pi*(wo-wt)/2);
        b_ap = [-eta -1];
        a_ap = [1 eta];
        [A,B,C,D] = tf2ss(b_ap,a_ap);
        if verbose
            disp('Appling balance transformation to the all-pass filter 1/F(z)');
        end
        [alpha,beta,gamma,delta] = balancing_transformation(A,B,C,D);
    case 'band'
        wb = omega(2);
        eta = cos(pi*wb);
        b_ap = [0 +eta -1];
        a_ap = [1 -eta 0];
        [A,B,C,D] = tf2ss(b_ap,a_ap);
        if verbose
            disp('Appling balance transformation to the all-pass filter 1/F(z)');
        end
        [alpha,beta,gamma,delta] = balancing_transformation(A,B,C,D);
end

% Obtain the state-space representation of the transformed filter
if verbose
    txt = ['Creating a composite maxflat ', type, '-pass filter G(z)'];
    disp(txt);
end
M = size(A_bal,1);
A_tf = kron(eye(M),alpha) + kron(A_bal*inv(eye(M)-delta.*A_bal),beta*gamma);
B_tf = kron(inv(eye(M)-delta.*A_bal)*B_bal,beta);
C_tf = kron(C_bal*inv(eye(M)-delta.*A_bal),gamma);
D_tf = D_bal + delta.*C_bal*inv(eye(M)-delta.*A_bal)*B_bal;

end

%==========================================================================
% function H = generate_filt_matrix_zero_inital(A,B,C,D,N)
%==========================================================================
% @description: Generates a lower-triangular Hankel matrix consisting of
%               the impulse response to the digital filter
% @inputs     :
%             - A, B, C, D - the state-space representation of the
%             transfomed filter
%             - N1 - The size of the input vector
%             - N2 - The size of K-sparse derivative (N-K)
% @outputs    :
%             - H - matrix consisting of the impulse response of the
%             composite filter
%==========================================================================        
function H = generate_filt_matrix_zero_inital(A,B,C,D,N1,N2)

h = zeros(N1,1);
A_old = A;
for i = 1:N1
    if i == 1
        h(i) = D;
    elseif i == 2
        h(i) = C*B;
    else
        h(i) = C*A_old*B;
        A_old = A_old*A;
    end   
end

% Use circulant shift property to compute the lower Hankel matrix
H = zeros(N1,N2);
for i = 1:N2
    H(:,i) = h;
    h = circshift(h,[1 0]);
    h(1:i) = 0;
end


end