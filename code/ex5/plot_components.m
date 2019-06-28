%% function plot_components( y, x_r, f_r, s_r, params, method )
%
% Inputs:
%   y       - noisy signal
%   x_r     - sparse and sparse derivative signal
%   f_r     - low-frequency signal
%   s_r     - oscillatory signal
%   params  - parameters for the labels of the plot
%   method  - SASDPR/DETOKS
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%

function plot_components( y, x_r, f_r, s_r, params, method )

txt_1 = [method, ', Noisy signal ($\mathbf{y}$)', ...
    ', $\lambda_1$ = ', num2str(params(1)), ', $\lambda_2$ = ', num2str(params(2)), ...
    ', $\lambda_3$ = ', num2str(params(3)), ', $\mu$ = ', num2str(params(4))];
txt_2 = [method, ', Low-frequency signal ($\mathbf{x}_1$)'];
txt_3 = [method, ', Sparse and sparse-derivative signal ($\mathbf{x}_2$)'];
txt_4 = [method, ', Oscillatory signal ($\mathbf{x}_3$)'];

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
plot(n,f_r,'k'); hold off;
title(txt_2,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,3);
text(-0.1,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,x_r,'k'); hold off;
title(txt_3,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,4);
text(-0.1,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,s_r,'k'); hold off;
title(txt_4,'interpreter','latex')
set(gca, 'box', 'off')
xlabel('Time (samples)','interpreter','latex')
axis tight;

printme_pdf = @(ex,meth) print('-dpdf', sprintf('../../results/%s_%s',ex,meth));
printme_pdf('ex5_1_',lower(method));

end

