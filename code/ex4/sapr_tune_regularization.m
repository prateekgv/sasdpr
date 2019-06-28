function [lam1, lam2, lam3] = sapr_tune_regularization(y, f, s, x, fs, fc, Nit, mu)

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
lam2 = 0.1:0.1:0.5;
lam3 = 0.1:0.1:0.5;
r = 1; P = fs/5; N = length(y);
y1 = preproc(r, P, y);

for ii = 1:length(lam1)
    for jj = 1:length(lam2)
        for kk = 1:length(lam3)
            [x_sapr,s_sapr,f_sapr,cost_sapr] = SAPR_v1(y1,fs,fc,lam1(ii),lam2(jj),lam3(kk),Nit,mu);
            x_sapr = x_sapr(P+1:P+N);
            s_sapr = s_sapr(P+1:P+N);
            f_sapr = f_sapr(P+1:P+N);
            rmse_x = sqrt(mean((x - x_sapr').^2));                        % rmse_x : rmse of sparse signal
            rmse_f = sqrt(mean((f - f_sapr').^2));                        % rmse_f : rmse of low-freqnecy signal
            % fprintf('lam1 = %.3f \t lam2 = %.3f \t lam3 = %.3f \t rmse_x = %.3f \t rmse_f = %.3f \n', lam1(ii),lam2(jj),lam3(kk),rmse_x,rmse_f);
            fprintf('lam1 = %.3f \t lam2 = %.3f \t lam3 = %.3f \t rmse = %.3f \n', lam1(ii),lam2(jj),lam3(kk),rmse_x+rmse_f);
        end
    end
end

end

