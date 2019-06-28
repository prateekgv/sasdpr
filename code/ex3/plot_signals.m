%% function plot_signals( s, y, x1, x2, x, method, params )
%
% Inputs:
%   s       - original signal
%   y       - noisy signal
%   x1      - reconstructed low-frequency signal
%   x2      - reconstructed k-order sparse derivative
%   x       - reconstructed signal
%   method  - SASS/SASD
%   params  - regularization and other params
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%

function plot_signals( y, x1, x2, x, method, params )

N = length(y);
n = 0:N-1;

txt_1 = ['Real ECG data'] ;
%txt_3 = sprintf('Sparse derivative (x2)');
txt_3 = ['Sparse derivative ($\mathbf{x}_2 = \mathbf{S} \mathbf{v}$)'];
%txt_2 = sprintf('Low-pass frequency (x1)');
txt_2 = ['Low-frequency signal ($\mathbf{x}_1$)'];
%txt_4 = sprintf('SAPR (K = %d, lambda = %.2f, M = %d, wc = %.3f),  RMSE = %.3f', K, lam, deg(1), fc, RMSE);
txt_4 = [method, ' (K = ', num2str(params(1)), ', $\lambda$ = ', num2str(params(2)), ', M = ', ...
    num2str(params(3)), ', $\omega_0$ = ',num2str(params(4)), '), RMSE = ',num2str(params(5),'%.3f')];
% txt_5 = ['Validation of optimality ($\ell_1$ norm penalty)'];

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
plot(n,x1,'k')
title(txt_2,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,3)
text(-0.1,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,x2,'k')
title(txt_3,'interpreter','latex')
set(gca, 'box', 'off')
axis tight;

subplot(4,1,4);
text(-0.1,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n,y,'k-.'); hold on;
plot(n,x,'k'); 
legend('Input','Reconstructed','location','north')
legend boxoff
title(txt_4,'interpreter','latex')
xlabel('Time (samples)','interpreter','latex')
set(gca, 'box', 'off')
axis tight;

printme_pdf = @(ex,meth) print('-dpdf', sprintf('../../results/%s_%s_%d',ex,meth));
printme_pdf('ex3',lower(method));

end

