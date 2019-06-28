%% Example: SASDPR (Spindle Detection)
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

%% Add paths
addpath('../detoks/');
addpath('../sapr/');
addpath('../utils/');

%% Load signal
load('../../data/spindles/excerpt5_ep28_spindle.mat');                         % load from excerpt6.edf epoch #33
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
tk_detoks = teager_operator(s_detoks);                                      % TKEO operator
[~,~,binDetoks] = detect_roi( tk_detoks, fs, th, 'spindle', 'detoks');
plot_signals(y, x_ex3, s_detoks, tk_detoks, th, binExpert, binDetoks, 'DETOKS');
plot_components(y, x_detoks, f_detoks, s_detoks, [lam0,lam1,lam2,mu], 'DETOKS');

%% Apply SAPR algorithm
% Parameters
lam0 = simdata.lam0_sapr;
lam1 = simdata.lam1_sapr;
lam2 = simdata.lam2_sapr;
nit = simdata.nit_sapr;
mu = simdata.mu_sapr;
th = simdata.th_sapr;
Hz = simdata.Hz;
fc = Hz/(fs/2);
P = fs/5; r = 1;

y1 = preproc(r, P, y);
[x_sapr,s_sapr,f_sapr,~,cost_sapr] = SASDPR_v1(y1,fs,fc,lam0,lam1,lam2,nit,mu);
x_sapr = x_sapr(P+1:P+N);
s_sapr = s_sapr(P+1:P+N);
f_sapr = f_sapr(P+1:P+N);
tk_sapr = teager_operator(s_sapr);
[~,~,binSapr] = detect_roi( tk_sapr, fs, th, 'spindle', 'sapr');
plot_signals(y, x_ex3, s_sapr, tk_sapr, th, binExpert, binSapr, 'SASDPR');
plot_components(y, x_sapr, f_sapr, s_sapr, [lam0,lam1,lam2,mu], 'SASDPR');