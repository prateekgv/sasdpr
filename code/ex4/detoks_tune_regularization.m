function [lam1, lam2, lam3] = detoks_tune_regularization(y, f, s, x, fs, deg, fc, Nit, mu)

% y - noisy signal
% f - low-frequency signal
% s - oscillatory signal
% x - sparse and sparse-derivative signal
% fs - sampling rate
% deg - deg of the filter
% fc - low-frequency cutoff signal
% Nit - no. of iterations
% mu - convergence rate

lam1 = 0.01:0.02:0.1;
lam2 = 0.1:0.2:1.0;
lam3 = 0.1:0.2:1.0;

for ii = 1:length(lam1)
    for jj = 1:length(lam2)
        for kk = 1:length(lam3)
            [x_detoks,s_detoks,f_detoks,cost_detoks] = DETOKS(y,fs,deg,fc,lam1(ii),lam2(jj),lam3(kk),Nit,mu);
            rmse_x = sqrt(mean((x - x_detoks').^2));                        % rmse_x : rmse of sparse signal
            rmse_f = sqrt(mean((f - f_detoks').^2));                        % rmse_f : rmse of low-freqnecy signal
            fprintf('lam1 = %.2f \t lam2 = %.2f \t lam3 = %.2f \t rmse_x = %.2f \t rmse_f = %.2f \n', lam1(ii),lam2(jj),lam3(kk),rmse_x,rmse_f);
        end
    end
end

end