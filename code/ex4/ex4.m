%% Example: SASDPR
% The sparsity-assisted signal denoising and pattern recognition (SASDPR)
% will simultaneously denoise and detect oscillatory pattern of interest in
% the given signal.
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019


%% Start
clear
close all
clc

%% Save fig command
printme = @(filename) print('-dpdf', sprintf('figures/Example4_%s', filename));

%% Add path
addpath('../sass/');
addpath('../detoks/');
addpath('../sapr/');
addpath('../utils/');

%% Create synthetic signal
% f - low-pass signal
% s - band-pass oscillating sparse signal
% x - high-pass sparse signal with sparse derivative
fs = 100; sigma = 0.1;
[y, f, s, x, ~] = generate_signal( fs, sigma );
N = length(y);
n = 0:N-1;

%% Apply DETOKS algorithm
% Parameters
Hz = 0.1;                                                                   % cutoff frequency of low-freq signal (Hz)
lam1_detoks = 0.05;                                                         % regularization on the sparsity of x
lam2_detoks = 0.50;                                                         % regularization on the sparsity of Dx
lam3_detoks = 0.15;                                                         % regularization on the sparsity of oscillatory components
d = 2;                                                                      % degree of the high-pass filter
fc = Hz/(fs/2);                                                             % normalized cutoff frequency of the low-freq signal
Nit = 50;                                                                   % number of iterations for the DETOKS algorithm
mu = 1.0;                                                                   % convergence rate parameter
[x_detoks,s_detoks,f_detoks,cost_detoks] = DETOKS(y,fs,d,fc,lam1_detoks,lam2_detoks,lam3_detoks,Nit,mu);
plot_signals( y, x, f, x_detoks, f_detoks, s_detoks, [lam1_detoks,lam2_detoks,lam3_detoks,mu,sigma],'DETOKS' )

%% SAPR parameters
lam1_sapr = 0.05;                                                           % regularization on the sparsity of x
lam2_sapr = 0.50;                                                           % regularization on the sparsity of Dx
lam3_sapr = 0.15;                                                           % regularization on the sparsity of oscillatory components
P = fs/5; r = 1;
y1 = preproc(r, P, y);
[x_sapr,s_sapr,f_sapr,cost_sapr] = SASDPR_v1(y1,fs,2*fc,lam1_sapr,lam2_sapr,lam3_sapr,Nit,mu);
x_sapr = x_sapr(P+1:P+N);
s_sapr = s_sapr(P+1:P+N);
f_sapr = f_sapr(P+1:P+N);
plot_signals( y, x, f, x_sapr, f_sapr, s_sapr,[lam1_sapr,lam2_sapr,lam3_sapr,mu,sigma],'SASDPR' )

