%% function plot_filter_response(H,b,a,N,omega,type,zero_phase)
%
% Inputs:
%   H       - Filter as matrix
%   b, a    - Filter coefficients
%   N       - Length of the filter
%   omega   - [wc, wn] - transformed cutoff frequency, original cutoff freq
%   type    - low, high, band pass filter
%   zero_phase - if true, zero phase filter
%
% Outputs:
%   None
%
% Prateek Gundannavar, prateek@wustl.edu, 2018
% Revised 2019
%% ________________________________________________________________________
%%


function plot_filter_response(H,b,a,N,omega,type,zero_phase)

%% save fig command
printme_type = @(zp,type,filename) print('-dpdf', sprintf('../../results/demo_%sbal_%s_%s',zp,type,filename));
% printme_eps = @(zp,type,filename) print('-depsc', sprintf('figures/demo_%sbal_%s_%s',zp,type,filename));

%% check if zero phase
if zero_phase
    mag = 1/2;
    H = H'*H;
    zp = 'zero_phase_';
else
    mag = 1/sqrt(2);
    zp = [];
end

%% Impulse input
x = zeros(N,1);
spike = N/2;
x(spike) = 1;

figure; clf;
subplot(2,2,1)
[Hf, om] = freqz(b,a);       % Fourier spectrum
if (zero_phase)
    plot(om/pi, abs(Hf).^2,'k','LineWidth',0.5);
else
    plot(om/pi, abs(Hf),'k','LineWidth',0.5);
end

if (strcmp(type,'high') || strcmp(type,'low'))
    line(omega(2), mag, 'marker', 'o','MarkerEdgeColor','black')
elseif (strcmp(type,'band'))
    line(omega(2)-omega(1)/2, mag, 'marker', 'o','MarkerEdgeColor','black'); hold on;
    line(omega(2)+omega(1)/2, mag, 'marker', 'o','MarkerEdgeColor','black'); hold off;
end

title('Frequency response')
legend(sprintf('%s-pass filter',type),'location','southeast')
xlabel('Normalized frequency')
ylim([0 1.2]); xlim([0,1.05]);
set(gca, 'box', 'off')
set(gca,'xtick',0:0.25:1); 
set(gca,'ticklabelinterpreter', 'tex');
set(gca,'xticklabel',{'0','0.25\pi','0.50\pi','0.75\pi','1.0\pi'});
set(gca,'ytick',0:0.5:1); 
ax = gca;
ax.XRuler.Axle.LineWidth = 0.5;
ax.YRuler.Axle.LineWidth = 0.5;

subplot(2,2,2)

if (zero_phase)
    b = conv(b,fliplr(b));
    a = conv(a,fliplr(a));
    [hz, hp, hl] = zplane(b, a);
else
    [hz, hp, hl] = zplane(b, a);
end
set(findobj(hz, 'Type', 'line'), 'Color', 'k'); 
set(findobj(hp, 'Type', 'line'), 'Color', 'k');
set(findobj(hl, 'Type', 'line'), 'Color', 'k');
set(gca, 'box', 'off')
title('Pole-zero diagram')

subplot(2,2,[3,4])
stem(1:N,H*x,...
    'Color','black',...
    'MarkerFaceColor','black',...
    'MarkerSize',3,...
    'MarkerEdgeColor','black');
xlabel('Time (samples)');
set(gca, 'box', 'off')
title(sprintf('Impulse response (n = %d)',spike));

%% Save plots
% printme_eps(zp,type,'filt_zplot')
printme_type(zp,type,'filt_zplot')
end
