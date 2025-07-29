% start with T (50001 points) and T2000 (2501 points)
% either loaded from documents, or calculated just now
% No equilibrium temperature calculated at this time
clear;
load('newresult2.mat');
%load('obsTemp2000.mat');
%T2000=Tobs;
T2000=T2000-mean(T2000);

dt=2000;    % 2000 years per point
Fs=1/dt;    % Sampling frequency 
N=length(T2000);    % N=2501
%N=1000;
frequencies=(0:N-1)*Fs/N;   % spectrum frequency
t=0:dt:5000000;
t1=0:100:5000000;
%===========================================================
% time series and spectrum

% run_ave=movmean(T,200);     % 20k year average
% run_jump=movmean(T,5000);   % 500k year average
run_ave=movmean(T2000,10);     % 20k year average using T2000
run_jump=movmean(T2000,250);    % 500k year average using T2000

%-----------------------------------------------------------
% Compute the FFT
fft_result=fft(T2000);
P=abs(fft_result).^2;
P=P/sum(P);

P2 = abs(fft_result/N).^2;         % two-sided power
P1 = P2(1:floor(N/2)+1);
P1(2:end-1) = 2*P1(2:end-1)/Fs*N;  % single-sided power
%-----------------------------------------------------------
% Welch’s power spectral density estimate
window=hamming(256);    % Hamming window of length 256
noverlap=128;   % 50% overlap
nfft=1024;   % Number of FFT points, can also be 256, 1024, etc.
[pxx,fwelch]= pwelch(T2000,window,noverlap,nfft,Fs);

%-----------------------------------------------------------
%figure;
% plot(t,T2000,'c',...
%     t,run_ave,t,run_jump);
% xlabel('Time (years)');
% ylabel('temperature (K)');
% legend({'T','movemean 1','movemean 2'},'Location','southwest');
% grid on;

% subplot(2,1,1)
% plot(t,T2000, 'LineWidth',0.2,'Color',[0.5 0.5 0.5])
% hold on
% plot(t,run_ave,'LineWidth',1.25)
% hold on
% plot(t,run_jump,'LineWidth',3,'Color',[1 0.8 0 0.5])
% hold off
% %xlabel('Time (10^6 years ago)');
% %ylabel('temperature (K)');
% legend({'T','20 kyr movemean','500 kyr movemean'},'Location','southwest');
% ylim([281,291.5])
% 
% ax=gca;
% ax.XTickLabel=flip(0:1:5);
% %ax.XAxis.Exponent=6;
% ax.FontSize = 14;
% 
% grid on;
% subplot(3,1,2)
% %loglog(frequencies, frequencies.*abs(fft_result));
loglog(frequencies(1:1250), P1(1:1250),'LineWidth',0.2,'Color',[0.5 0.5 0.5]);
%title('Spectrum');
%xlabel('Frequency (year^{-1})');
%ylabel('normalized power');
xlim([0,0.5*Fs]);
grid on;
hold on;
% subplot(3,1,3)
% %plot(fwelch, 10*log10(pxx).*fwelch);
loglog(fwelch,pxx,'LineWidth',1.25);
% xlabel('Frequency (year^-1)');
% ylabel('psd (K^2 year)');
% xlim([0,max(fwelch)])
% grid on;
hold off
ax=gca;
ax.FontSize = 14;
legend({'Fourier','Welch'},'Location','southwest');
%============================================================
% telegraph approximation

% teleAprox(T2000>run_jump)=1;
% interval=find(diff(teleAprox) ~= 0);
% interval(2:end)=diff(interval);
% tau_high=interval(1:2:end); % Elements at odd indices
% if (mod(length(interval),2)==1)
%     tau_high=tau_high(1:end-1);
% end
% tau_low=interval(2:2:end); % Elements at even indices
% tau_period=tau_high+tau_low;
% binEdge=min(tau_high):max(tau_high);
% [counts,edges]=histcounts(tau_high,binEdge, 'Normalization', 'pdf');
% %------------------------------------------------------------
% figure;
% subplot(2,1,1)
% area(t,teleAprox,'LineStyle','none');
% xlabel('Time (years)');
% ylim([0,1.2]);
% title('telegraph approximation');
% subplot(2,1,2)
% scatter(edges(1:end-1)*dt,counts);
% xlabel('tau(years)');
% %xlim([0,5E4])
% title('Probability Density');
% set(gca, 'YScale', 'log');
% 
% %=============================================================
% % auto-corr
% trans1=1340;     %find(run_ave<(0.5*T_unstable+0.5*T_high),1);
% trans2=1980;     %find(run_ave<T_unstable,1);
% [autoc0,lag0]=autocorr(T2000,N-1);
% [autoc1,lag1]=autocorr(T2000(1:trans1),trans1-1);
% [autoc2,lag2]=autocorr(T2000(trans1:trans2),trans2-trans1);
% [autoc3,lag3]=autocorr(T2000(trans2:end),N-trans2);
% tauId0=find(autoc0<0,1);
% tau0=sum(autoc0(1:tauId0-1))*dt;
% tauId1=find(autoc1<0,1);
% tau1=sum(autoc1(1:tauId1-1))*dt;
% tauId2=find(autoc2<0,1);
% tau2=sum(autoc2(1:tauId2-1))*dt;
% tauId3=find(autoc3<0,1);
% tau3=sum(autoc3(1:tauId3-1))*dt;
% %--------------------------------------------------------------
% figure;
% subplot(2,3,[1,3])
% plot(t,autoc0,t(1:trans1),autoc1,t(1:trans2-trans1+1),autoc2,t(1:end-trans2+1),autoc3);
% ylim([-0.5,1]);
% title('autocorrelation')
% legend('whole series','period 1','period 2','period 3')
% hold on
% xline(tau0);
% xlabel({'tau=',tau0});
% hold off
% grid on;
% 
% subplot(2,3,4)
% plot(t(1:trans1),autoc1);
% xlim([0,trans1*dt]);
% ylim([-0.5,1]);
% hold on
% xline(tau1);
% xlabel({'tau=',tau1});
% hold off
% grid on;
% 
% subplot(2,3,5)
% plot(t(1:trans2-trans1+1),autoc2);
% xlim([0,(trans2-trans1+1)*dt]);
% ylim([-0.5,1]);
% hold on
% xline(tau2);
% xlabel({'tau=',tau2});
% hold off
% grid on;
% 
% subplot(2,3,6)
% plot(t(1:end-trans2+1),autoc3);
% xlim([0,(N-trans2+1)*dt]);
% ylim([-0.5,1]);
% hold on
% xline(tau3);
% xlabel({'tau=',tau3});
% hold off
% grid on;

%=============================================================
% wavelet transform
[cwtcoeffs, wfreq] = cwt(T2000,'morse');
cwtf=sum(abs(cwtcoeffs).^2,2)/sum(abs(cwtcoeffs).^2,'all');

%-------------------------------------------------------------
%figure;
%subplot(1,3,[1,2])
% subplot(2,1,2)
% contourf(t, wfreq, abs(cwtcoeffs),0:0.05:2,'LineColor','none');
% ax=gca;
% ax.XTickLabel=flip(0:1:5);
% ax.XAxis.Exponent=6;
% ax.FontSize = 14;
% 
% 
% set(gca, 'YScale', 'log');
% 
% %period_ticks=[20.5,41,50,82,100,123,164,200];
% period_ticks=[800,400,200,100,41,20.5];
% freq_ticks=2./period_ticks;
% ax.YTick=freq_ticks;
% 
% 
% ax.YTickLabel = string(period_ticks);

% Define the anchor colors as RGB triplets
anchorColors = [
    1.00, 1.00, 1.00;   % white

    % Extended grey range
    0.98, 0.98, 0.98;
    0.95, 0.95, 0.95;
    0.91, 0.91, 0.91;
    0.87, 0.87, 0.87;
    0.82, 0.82, 0.82;
    0.79, 0.79, 0.79;
    0.82, 0.82, 0.79;
    0.87, 0.87, 0.79;
    0.93, 0.93, 0.79;


    % Extended light yellow range
    0.97, 0.98, 0.79;
    0.97, 0.98, 0.75;
    0.97, 0.98, 0.70;
    1.00, 0.96, 0.70;
    1.00, 0.96, 0.65;
    1.00, 0.93, 0.65;
    1.00, 0.93, 0.60;
    1.00, 0.93, 0.55;
    1.00, 0.90, 0.55;
    1.00, 0.90, 0.50;

    % Yellow to orange
    1.00, 0.85, 0.35;
    0.98, 0.75, 0.20;
    0.95, 0.65, 0.12;
    0.90, 0.55, 0.08;
    0.85, 0.45, 0.05;

    % Orange to bright red
    0.80, 0.36, 0.05;
    0.77, 0.30, 0.08;
    0.75, 0.25, 0.12;
    0.80, 0.22, 0.15;
    0.85, 0.20, 0.18;

    % Bright red shades
    0.90, 0.24, 0.22;
    0.95, 0.28, 0.27;
    0.98, 0.32, 0.31;
    1.00, 0.35, 0.35;

    % Pinkish reds to pastel pinks
    1.00, 0.50, 0.55;
    1.00, 0.57, 0.63;
    1.00, 0.65, 0.70;
    1.00, 0.75, 0.80;
    1.00, 0.83, 0.88;
    %1.00, 0.90, 0.93;

    
];

% % Number of colors desired in final colormap
% numColors = 40;
% 
% % Positions of anchor colors (scaled 0 to 1)
% x = linspace(0, 1, size(anchorColors,1));
% 
% % Positions for interpolation
% xi = linspace(0, 1, numColors);
% 
% % Interpolate each RGB channel
% R = interp1(x, anchorColors(:,1), xi);
% G = interp1(x, anchorColors(:,2), xi);
% B = interp1(x, anchorColors(:,3), xi);
% 
% % Combine into colormap matrix
% customColormap = [R' G' B'];

% Apply the colormap and show colorbar
%colormap(anchorColors);
%colorbar;
%xlabel('Time (10^6 years ago)');
%ylabel('period (kyr)');
%title('Magnitude Scalogram (CWT)');

% subplot(1,3,3)
% plot(wfreq,cwtf)
% set(gca, 'XScale', 'log');
% xlabel('Frequency(period=1/f*dt)')
% ylabel('variance')

% figure;
% plot(0:2000:5000000,abs(cwtcoeffs(33,:)),'LineWidth',1.5,'Color','#28557C')
% hold on
% plot(0:2000:5000000,abs(cwtcoeffs(45,:)),'LineWidth',1.5,'Color','#FA8072')
% hold off
% ylim([0,2.1])
% ax=gca;
% ax.XTickLabel=flip(0:1:5);
% ax.XAxis.Exponent=6;
% ax.FontSize = 14;
% %xlabel('Time (10^6 years ago)');
% %ylabel('amplitude');
% legend(["41 kyr","100 kyr"],'Location','northwest')
