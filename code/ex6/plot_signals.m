%% plot_signals( y, x_ex3, kc, k, tk, th, binX, binA, method)
%
% Inputs:
%   y       - noisy signal
%   x_ex3   - annotated sleep spindle signal
%   k       - K-complex signal detected
%   tk      - Teager-Kaiser energy operator
%   th      - threshold for spindle detection
%   binX    - binary vector with expert annotated spindle
%   binA    - binary vector with algorithm annotated spindle
%   method  - SASDPR/DETOKS
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%
function plot_signals( y, x_ex3, kc, k, tk, th, binX, binA, method)

global simdata;
fs = simdata.fs;

if strcmp(method, 'DETOKS')
    txt_1 = [method, ', EEG signal ($\mathbf{y}$), ','$\lambda_1$ = ', num2str(simdata.lam0_detoks), ...
        ', $\lambda_2$ = ', num2str(simdata.lam1_detoks), ', $\lambda_3$ = ', num2str(simdata.lam2_detoks),...
        ', $\mu$ = ', num2str(simdata.mu_detoks)];
    txt_2 = [method, ', Low-frequency signal'];
    txt_3 = [method, ', Teager-Kaiser energy operator, Threshold = ', num2str(simdata.th_detoks)];
    txt_4 = [method, ', Detected EEG K-complexes'];
elseif strcmp(method, 'SAPR')
    txt_1 = [method, ', EEG signal ($\mathbf{y}$), ','$\lambda_1$ = ', num2str(simdata.lam0_sapr), ...
        ', $\lambda_2$ = ', num2str(simdata.lam1_sapr), ', $\eta$ = ', num2str(simdata.eta_sapr),...
        ', $\mu$ = ', num2str(simdata.mu_sapr)];
    txt_2 = [method, ', K-Complex signal ($\mathbf{x}_2$)'];
    txt_3 = [method, ', Teager-Kaiser energy operator, Threshold = ', num2str(simdata.th_sapr)];
    txt_4 = [method, ', Detected EEG K-complexes'];
end

N = length(y);
n = 0:N-1;

%% Plot
figure('rend','painters','pos',[100 100 550 500]);
clf

subplot(4,1,1);
text(-0.1,0.98,'(a)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n/fs, y, 'k'); hold on;
p0 = plot(n/fs, x_ex3, 'r'); hold off;
title(txt_1,'interpreter','latex')
legend([p0], {'K-complex'},'location','north');
legend boxoff;
set(gca, 'box', 'off')
axis tight;

subplot(4,1,2);
text(-0.1,0.98,'(b)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
p0 = plot(n/fs, kc, 'k'); hold on;
p1 = plot(n/fs, k, '--k'); hold off;
title(txt_2,'interpreter','latex')
if strcmp(method, 'SAPR')
    legend([p1,p0], {'$\mathbf{\Psi} \mathbf{k}$', ...
        '$\mathbf{B}^T \mathbf{B} \mathbf{\Psi} \mathbf{k}$'},'interpreter','latex', ...
        'location','north');
    legend boxoff;
end
set(gca, 'box', 'off')
axis tight;

subplot(4,1,3);
text(-0.1,0.98,'(c)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
p0 = plot(n/fs, tk); hold on;
p1 = plot(n/fs, th*ones(1,N), 'r--'); hold off;
title(txt_3,'interpreter','latex')
legend([p0,p1], {'TKEO','Threshold'},'location','north');
legend boxoff;
set(gca, 'box', 'off')
axis tight;

subplot(4,1,4);
text(-0.1,0.98,'(d)','Units', 'Normalized', 'VerticalAlignment', 'Top'); hold on;
plot(n/fs, binX, 'Color', [0,0,0]+0.7,'linewidth',1.0); hold on;
plot(n/fs, binA, 'k');
title(txt_4,'interpreter','latex')
xlabel('Time (seconds)','interpreter','latex')
legend('Expert',method,'location','north');
legend boxoff
set(gca, 'box', 'off')
ylim([-0.2, 1.2])

printme_pdf = @(ex,meth) print('-dpdf', sprintf('../../results/%s_%s',ex,meth));
printme_pdf('ex6',lower(method));

% printme_eps = @(ex,meth) print('-depsc', sprintf('figures/%s_%s',ex,meth));
% printme_eps('ex6',lower(method));


end

