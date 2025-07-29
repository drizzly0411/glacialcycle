clear;

C=4.2E8;
tspan = [0 5000000]; 
T0 = 289;       

dt=100;
t = tspan(1):dt:tspan(2);
T = zeros(size(t));
T_unstable=T;
T_high=T;
T_low=T;
teleAprox=T;
T(1) = T0;

noiseWhite=dsp.ColoredNoise('Color','white', 'SamplesPerFrame',length(t)-1);
noiseBlue=dsp.ColoredNoise('Color','blue', 'SamplesPerFrame',length(t)-1);
noisePink=dsp.ColoredNoise('Color','pink', 'SamplesPerFrame',length(t)-1);
noise=8*noiseWhite()+0*noiseBlue()+0*noisePink();
std0=std(noise);
noise=(noise-mean(noise))/std0;

Am=0.05;  % amplitude of the noise
for i = 2:length(t)
    [dT_dt0,high,unstable,low] = dTdt3(t(i-1), T(i-1));
    % if t(i-1)>2000000
    %     Am=Am-0.02/30000;%0.09-0.04*(t(i-1)/2500000-1);
    % end
    dT_dt=dT_dt0+Am*noise(i-1)/C;
    T(i)=T(i-1)+dt*dT_dt*86400*365;
    T_unstable(i)=unstable;
    T_high(i)=high;
    T_low(i)=low;
    %noise(i-1)=noise(i-1)*Am;
end
%T(2:end)=T0+noise*Am;
%plot(noise)


%------ large ensemble ------

% NEns=1000;
% TEns=T/NEns;
% for jj=2:NEns
%     noise=8*noiseWhite()+0*noiseBlue()+0*noisePink();
%     std0=std(noise);
%     noise=(noise-mean(noise))/std0;
%     Am=0.05;
%     for i = 2:length(t)
%         % if t(i-1)>2000000
%         %     Am=Am-0.02/30000; %0.09-0.04*(t(i-1)/2500000-1);
%         % end
%         [dT_dt0,~,~,~] = dTdt3(t(i-1), T(i-1));
%         dT_dt=dT_dt0+Am*noise(i-1)/C;
%         T(i)=T(i-1)+dt*dT_dt*86400*365;
%     end
% 
%     TEns=TEns+T/NEns;
% 
% end
% T=TEns;
%----------------------------

run_ave=movmean(T,200);
run_jump=movmean(T,5000);

T_unstable(1)=T_unstable(2);
T_high(1)=T_high(2);
T_low(1)=T_low(2);

%==============
T2000=zeros([1,2501],'double');
rshp=reshape(T(2:end),20,[]);
T2000(2:end)=mean(rshp);
T2000(1)=T(1);

figure;
%scatter(t,T_high,0.1)
%hold on
%scatter(t,T_low,0.1)
%hold on
%scatter(t,T_unstable,0.1)
%hold on


% plot(tspan(1):dt:tspan(2),T)
% hold on
%plot(tspan(1):dt*20:tspan(2),T2000)
%hold on
%plot(t,run_jump,t,run_ave)
%hold off
%xlim([0,50000])

[cwtcoeffs, wfreq] = cwt(T2000,'morse');
cwtf=sum(abs(cwtcoeffs).^2,2)/sum(abs(cwtcoeffs).^2,'all');

%-------------------------------------------------------------

subplot(3,1,1)
plot(tspan(1):dt*20:tspan(2),T2000)
% hold on
% plot(t,run_jump,t,run_ave)
hold off

subplot(3,1,[2,3])
contourf(0:2000:5000000, wfreq, abs(cwtcoeffs),0:0.05:1.8,'LineColor','none');
%contourf(0:2000:5000000, wfreq, abs(cwtcoeffs),'LineColor','none');
set(gca, 'YScale', 'log');
colorbar;
xlabel('Time (years)');
ylabel('Frequency');
title('Magnitude Scalogram (CWT)');

% subplot(2,3,6)
% plot(wfreq,cwtf)
% set(gca, 'XScale', 'log');
% xlabel('Frequency(period=1/f*dt)')
% ylabel('variance')

% figure;
% plot(0:2000:5000000,abs(cwtcoeffs(33,:)))
% hold on
% plot(0:2000:5000000,abs(cwtcoeffs(45,:)))
% hold off
% %ylim([0,1.8])
% ax=gca;
% ax.XTickLabel=flip(0:0.5:5);
% ax.XAxis.Exponent=6;
% xlabel('Time (10^6 years ago)');
% ylabel('amplitude');
% legend(["41 kyr","100 kyr"])

save('temporary.mat','T','T2000');

