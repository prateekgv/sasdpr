%% Example: SASD
% The sparsity-assisted signal denoising (SASD) algorithm simultaneously
% combines total variation denoising and low-pass filtering to filter a
% noisy signal, such that it preserves the discontinuities. In this
% example, we denoise ECG signal using LPF, TVD, SASS, and SASD
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019

%% Start
clear
close all
clc

%% Add path
addpath('../sass/');
addpath('../sasd/');
addpath('../utils/');

%% Load precomputed factorized matrices for SASD
preload = true;                                                             % load precomputed matrices
show_plots = true;                                                          % display plots

%% Make signal
load('../../data/108m_2_3.mat');                                               % load 'val'

N = 1000;                                                                   % length of signal to extract
dat = val(1, 2360+(1:N));                                                   % extract signal
dat = dat(:);                                                               % convert to column vector
y = dat - mean(dat);                                                        % remove mean

n = 1:N;
figure('rend','painters','pos',[100 100 550 500]);
clf

%% Low-pass filtering
fc = 0.02;                                                                  % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
d = 3;                                                                      % d : filter order parameter (d = 1, 2, or 3)
K = 2;                                                                      % K : order of sparse derivative
[A, B, ~, ~, ~, ~, ~, ~, ~] = ABfilt(d, fc, N, K);                          % Sparse-banded filter matrices
L = @(x) (x - A\(B*x));                                                     % L : low-pass filtering
x_lpf = L(y);                                                               % Apply low-pass filtering
rmse_lpf = sqrt(mean((y - x_lpf).^2));                                      % rmse_lpf : rmse of low-pass filtering

subplot(4,1,1);
text(-0.1,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'color', [0.7,0.7,0.7],'linewidth',1.0); hold on;
plot(n,x_lpf,'k'); 
legend('Input','Reconstructed','location','northwest')
legend boxoff
txt_1 = ['LPF', ' ($\omega_0$ = ', num2str(fc), ')'];
title(txt_1,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

%% TV denoising
Nit = 100;
lam_tvd = 20;
[x_tvd, cost] = tvd_mm(y, lam_tvd, Nit);
rmse_tvd = sqrt(mean((y - x_tvd).^2));                                      % rmse_lpf : rmse of low-pass filtering

subplot(4,1,2);
text(-0.1,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'color', [0.7,0.7,0.7],'linewidth',1.0); hold on;
plot(n,x_tvd,'k'); 
legend('Input','Reconstructed','location','northwest')
legend boxoff
txt_1 = ['TVD', ' ($\lambda$ = ', num2str(lam_tvd), ')'];
title(txt_1,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

%% SASS parameters
fc = 0.02;                                                                  % fc : cut-off frequency (cycles/sample) (0 < fc < 0.5);
d = 3;                                                                      % d : filter order parameter (d = 1, 2, or 3)
K = 2;                                                                      % K : order of sparse derivative
lam_sass = 200;                                                             % lam_sass : regularization parameter
beta = 1;                                                                   % scaling for regularization
[A, B, B1, D, ~, ~, ~, ~, HTH1norm_sass] = ABfilt(d, fc, N, K);             % Sparse-banded filter matrices
tic
[x_sass, ~] = sass_L1(y, d, fc, K, lam_sass);                               % Apply SASS
toc
rmse_sass = sqrt(mean((y - x_sass).^2));                                    % rmse_sass : rmse of SASS

subplot(4,1,3);
text(-0.1,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'color', [0.7,0.7,0.7],'linewidth',1.0); hold on;
plot(n,x_sass,'k'); 
legend('Input','Reconstructed','location','northwest')
legend boxoff
txt_3 = ['SASS', ' (K = ', num2str(K), ', $\lambda$ = ', ...
    num2str(lam_sass), ', M = ', num2str(2*d), ', $\omega_0$ = ', ...
    num2str(fc), ')'];
title(txt_3,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

%% SASD parameters
deg = [4,4];                                                                % deg : filter order parameter (numerator and denominator)
wc  = fc*2;                                                                 % wc : cut-off frequency (cycles/sample) (0 < wc < 1);
zp = 100;                                                                   % zp : length of preprocessing
r = 1;                                                                      % r : degree of the polynomial fit for the preprocessed signal
lam_sasd = 200;                                                             % lam_sasd : regularization parameter
wn = 0.1;                                                                   % fn : prototype low-pass filter cutoff frequency
if ~preload
    [H,H1,~,~,~,H2norm,~] = IIR_ABfilt(deg, N+2*zp, [wn,wc], 'high', K);    % Filters as matrices (not sparse)
    save('../../data/precomputed_mats_3.mat','H','H1','H2norm');
else
    disp('Loading precomputed matrices as filters');
    load '../../data/precomputed_mats_3.mat'
end

y1 = preproc(r, zp, y);                                                     % Preprocess the signal
tic
[x_sasd, ~, ~, ~, ~] = sasd_L1(y1, K, lam_sasd, H, H1);                     % Apply SASD
toc
x_sasd = x_sasd(zp+1:zp+N);                                                 % Remove the sides
rmse_sasd = sqrt(mean((y - x_sasd).^2));                                    % rmse_sasd : rmse of SASD

subplot(4,1,4);
text(-0.1,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'color', [0.7,0.7,0.7],'linewidth',1.0); hold on;
plot(n,x_sasd,'k'); 
legend('Input','Reconstructed','location','northwest')
legend boxoff
txt_4 = ['SASD', ' (K = ', num2str(K), ', $\lambda$ = ', ...
    num2str(lam_sasd), ', M = ', num2str(2*deg(1)), ', $\omega_0$ = ', ...
    num2str(fc), ')'];
title(txt_4,'interpreter','latex')
xlabel('Time (samples)','interpreter','latex')
set(gca, 'box', 'off')
axis tight;

%% RMSE
printme_pdf = @(ex,meth) print('-dpdf', sprintf('../../results/%s_%s',ex,meth));
printme_pdf('ex3','compare');