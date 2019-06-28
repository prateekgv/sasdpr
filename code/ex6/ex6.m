%% Example: SAPR (K-complex Detection)
% The sparsity-assisted pattern recognition (SAPR) will detect
% wavelet-like patterns in the signal.
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019

%% Start
clear
close all
clc

%% Add path
addpath('../detoks/');
addpath('../sapr/');
addpath('../utils/');

%% Make signal
load('../../data/kcomplexes/excerpt1_ep18_kcomplex.mat');                           % load from excerpt6.edf 
load_settings();

global simdata;
fs = simdata.fs;                                                            % sampling rate of EEG
y = x_in;
N = length(y);
n = 0:N-1;
binExpert = zeros(1,N);
binExpert(~isnan(x_ex3)) = 1;                                               % Expert detected spindles

%% Apply DETOKS algorithm
% Parameters
lam0 = simdata.lam0_detoks;
lam1 = simdata.lam1_detoks;
lam2 = simdata.lam2_detoks;
deg = simdata.deg;
nit = simdata.nit_detoks;
mu = simdata.mu_detoks;
th = simdata.th_detoks;
fc = 1/(fs/2);

[x_detoks,s_detoks,f_detoks,cost_detoks] = DETOKS(y,fs,deg,fc,lam0,lam1,lam2,nit,mu);
tk_detoks = teager_operator(f_detoks);                                      % TKEO operator
[~,~,binDetoks] = detect_roi( tk_detoks, fs, th, 'kcomplex', 'detoks');
plot_signals(y, x_ex3, f_detoks, f_detoks, tk_detoks, th, binExpert, binDetoks, 'DETOKS');

%% Apply SAPR algorithm
% Parameters
lam0 = simdata.lam0_sapr;
lam1 = simdata.lam1_sapr;
nit = simdata.nit_sapr;
mu = simdata.mu_sapr;
eta = simdata.eta_sapr;
th = simdata.th_sapr;
Hz = simdata.Hz;
fc = Hz/(fs/2);
P = fs/5; r = 1;

y1 = preproc(r, P, y);
[k_bp_sapr,k_sapr,cost_sapr] = SAPR_v1(y1,fs,fc,lam0,lam1,mu,eta,nit);
k_bp_sapr = k_bp_sapr(P+1:P+N);
k_sapr = k_sapr(P+1:P+N);
tk_sapr = teager_operator(k_bp_sapr);
[~,~,binSapr] = detect_roi( tk_sapr, fs, th, 'kcomplex', 'sapr');
plot_signals(y, x_ex3, k_bp_sapr, k_sapr, tk_sapr, th, binExpert, binSapr, 'SAPR');