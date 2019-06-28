%==========================================================================
% @mat-script   : demo_filt_matrices.m
% @author       : Prateek Gundannavar
% @description  : An example to demonstrate zero-phase filtering using
%                 banded matrices 
%==========================================================================

%%  Start
close all; clear all; clc;

%% Add path
addpath('../sasd/');

%% Define prototype low-pass filter using maxflat.m
NA = 2;                     % NA : The degree of the denominator polynomial
NB = 2;                     % NB : The degree of the numerator polynomial
deg = [NA, NB];             % deg : The degree of the numerator and denominator polynomial
wn = 0.1;                   % wn : The cut-off frequency of the prototype filter
N = 100;                    % N  : Size of the matrix
zero_phase = true;             % causal : Indicates if the response is causal or not
K = 1;  

%% Define the cut-off and central frequency of composite LPF
wc = 0.2;                   % wc : The cut-off frequency of the composite filter
type = 'low';
[L,L_lp,D_lp,b_lp,a_lp,~,~] = IIR_ABfilt(deg,N,[wn,wc],type);
plot_filter_response(L,b_lp,a_lp,N,[wn,wc],type,zero_phase)
plot_filter_response(L,b_lp,a_lp,N,[wn,wc],type,~zero_phase)

%% Define the cut-off and central frequency of composite HPF
wc = 0.2;                   % wc : The cut-off frequency of the composite filter
type = 'high';
[H,H_hp,D_hp,b_hp,a_hp,~,err] = IIR_ABfilt(deg,N,[wn,wc],type,K);
fprintf('Sparsity K = %d, Norm error = %.3f \n',K,err);
plot_filter_response(H,b_hp,a_hp,N,[wn,wc],type,zero_phase)
plot_filter_response(H,b_hp,a_hp,N,[wn,wc],type,~zero_phase)

%% Define the bandwidth and central frequency of composite BPF
wb = 0.5;                   % wb : The central frequency
type = 'band';
[B,B_bp,D_bp,b_bp,a_bp,~,err] = IIR_ABfilt(deg,N,[wn,wb],type,K);
fprintf('Sparsity K = %d, Norm error = %.3f \n',K,err);
plot_filter_response(B,b_bp,a_bp,N,[wn,wb],type,zero_phase)
plot_filter_response(B,b_bp,a_bp,N,[wn,wb],type,~zero_phase)
