%% Example: SASD
% The sparsity-assisted signal denoising (SASD) algorithm simultaneously
% combines total variation denoising and low-pass filtering to filter a
% noisy signal, such that it preserves the discontinuities.
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019

%% Start
clear
close all
clc

%% Save fig commands
printme = @(filename) print('-dpdf', sprintf('figures/Example1_%s', filename));
addpath('../sass/');
addpath('../sasd/');
addpath('../utils/');


%% Load precomputed factorized matrices for SASD
preload = true;                                                             % load precomputed matrices
show_plots = true;                                                          % display plots

%% Create test signal
N = 300;                                                                    % N : length of signal
n = (1:N)';

s1 = 2*(n < 0.3*N) + 1*(n > 0.6*N);                                         % s1 : step function
s2 = sin(0.021*pi*n);                                                       % s2 : smooth function
s = s1 + s2;                                                                % s : total signal

%% SASS parameters
fc = 0.022;                                                                 % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
d = 3;                                                                      % d : filter order parameter (d = 1, 2, or 3)
K = 1;                                                                      % K : order of sparse derivative
[A, B, B1, D, ~, ~, ~, ~, HTH1norm_sass] = ABfilt(d, fc, N, K);             % Sparse-banded filter matrices

%% SASD parameters
deg = [3,3];                                                                % d : filter order parameter (numerator and denominator)
wc  = 0.022*2;                                                              % wc : cut-off frequency (cycles/sample) (0 < wc < 1);
zp = 20;                                                                    % zp : length of preprocessing
r = 2;                                                                      % r : degree of the polynomial fit for the preprocessed signal
wn = 0.1;                                                                   % fn : prototype low-pass filter cutoff frequency
if ~preload
    [H, H1, ~, ~, ~, H2norm, ~] = IIR_ABfilt(deg, N+2*zp, [wn,wc],...
        'high', K);                                                         % Filters as matrices (not sparse)
    save('precomputed_mats_1.mat','H','H1','H2norm');                         % Save precomputed matrices for future use
else
    disp('Loading precomputed matrices as filters')
    load '../../data/precomputed_mats_1.mat'
end

%% Create noisy signal
disp('Running fixed random seed');
randn('seed',1);                                                            % Set state so example is reproducible
sigma = 0.1:0.1:0.5;                                                        % sigma : standard deviation of noise

for ii = 1:length(sigma)
    
    %% SASS
    noise = sigma(ii)*randn(N+2*zp, 1);                                     % noise : white zero-mean Gaussian
    y = s + noise(zp+1:N+zp,1);                                             % y : noisy data
    lam = determine_lambda_sass(A, B, N, K, noise(zp+1:N+zp));
    % [x_sass, ~] = sass_L1(y, d, fc, K, 2*HTH1norm_sass*sigma(ii));        % Apply SASS
    [x_sass, ~] = sass_L1(y, d, fc, K, lam);          % Apply SASS
    rmse_sass = sqrt(mean((s - x_sass).^2));                                % rmse_sass : rmse of SASS
    
    %% SASD
    y1 = preproc(r, zp, y);                                                 % Preprocess the signal
    lam = determine_lambda_sasd(H, N+2*zp, K, noise);
    % [x_sasd, x1, x2, ~, ~] = sasd_L1(y1, K, 2*H2norm*sigma(ii), H, H1);   % Apply SASD
    [x_sasd, x1, x2, ~, ~] = sasd_L1(y1, K, lam, H, H1);     % Apply SASD
    x_sasd = x_sasd(zp+1:zp+N);                                             % Remove the sides
    rmse_sasd = sqrt(mean((s - x_sasd).^2));                                % rmse_sasd : rmse of SASD
    
    fprintf('Sigma = %f \t RMSE SASS: %f \t RMSE SASD: %f \n', ...
        sigma(ii),rmse_sass,rmse_sasd);
    
    if show_plots
        params = zeros(1,6);
        params(1) = sigma(ii); params(2) = K; params(3) = lam;
        params(4) = deg(1)*2; params(5) = fc; params(6) = rmse_sasd;
        x1 = x1(zp+1:zp+N); x2 = x2(zp+1:zp+N);
        plot_signals( s, y, x1, x2, x_sasd, 'SAPR', params );
        close all;
    end
end

