%% function plot_signals( y, x, f, x_r, f_r, s_r, params, method )
%
% Inputs:
%   y       - noisy signal
%   x       - original sparse signal
%   f       - original low-frequency signal
%   x_r     - reconstructed sparse signal
%   f_r     - reconstructed low-frequency signal
%   s_r     - oscillatory pattern
%   params  - regularization and other params
%   method  - SASDPR/DETOKS
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%

function plot_signals( y, x, f, x_r, f_r, s_r, params, method )

txt_1 = [method, ', Noisy signal ($\mathbf{y}$), ' '$\sigma = $ ', num2str(params(5)), ...
    ', $\lambda_1$ = ', num2str(params(1)), ', $\lambda_2$ = ', num2str(params(2)), ...
    ', $\lambda_3$ = ', num2str(params(3)), ', $\mu$ = ', num2str(params(4))];
rmse = sqrt(mean((f - f_r').^2));
txt_2 = [method, ', Low-frequency signal ($\mathbf{x}_1$), RMSE = ', num2str(rmse,'%.3f')];
rmse = sqrt(mean((x - x_r').^2));
txt_4 = [method, ', Sparse and sparse-derivative signal ($\mathbf{x}_3$), RMSE = ', num2str(rmse,'%.3f')];
txt_3 = [method, ', Oscillatory signal ($\mathbf{x}_2$)'];

N = length(y);
n = 0:N-1;

%% Plot
figure('rend','painters','pos',[100 100 550 500]);
clf

subplot(4,1,1);
text(-0.1,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'k')
title(txt_1,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,2);
text(-0.1,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,f,'Color', [0,0,0]+0.7,'linewidth',1.0); hold on;
plot(n,f_r,'k'); hold off;
legend('Input','Reconstructed','location','southeast')
legend boxoff
title(txt_2,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,3);
text(-0.1,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,s_r,'k'); hold off;
title(txt_3,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,4);
text(-0.1,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,x,'Color', [0,0,0]+0.7,'linewidth',1.0); hold on;
plot(n,x_r,'k'); hold off;
xlabel('Time (samples)','interpreter','latex')
legend('Input','Reconstructed','location','east')
legend boxoff
title(txt_4,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

printme_pdf = @(ex,meth) print('-dpdf', sprintf('../../results/%s_%s',ex,meth));
printme_pdf('ex4',lower(method));

% printme_eps = @(ex,meth) print('-depsc', sprintf('figures/%s_%s',ex,meth));
% printme_eps('ex4',lower(method));

end

