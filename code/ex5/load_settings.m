%% function load_settings()
%
% Description:
%   Load settings of the SASDPR and DETOKS algorithms
%
% Inputs:
%   None
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%


function load_settings()

global simdata;

%% General settings
simdata.fs = 200;
simdata.epoch = 30;

%% DETOKS settings
simdata.method = 'DETOKS';
simdata.nit_detoks = 50;
simdata.mu_detoks = 0.5;
simdata.lam0_detoks = 0.60;
simdata.lam1_detoks = 7.00;
simdata.lam2_detoks = 7.00;
simdata.deg = 2;
simdata.Hz = 1.0;
simdata.th_detoks = 0.03;

%% SAPR settings
simdata.method = 'SAPR';
simdata.nit_sapr = 50;
simdata.mu_sapr = 0.1;
simdata.lam0_sapr = 0.5;
simdata.lam1_sapr = 4.8;
simdata.lam2_sapr = 5.6;
simdata.Hz = 2.0;
simdata.th_sapr = 0.05;

end

