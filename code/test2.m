x=250:0.2:300;  % temperature (K)
xN=length(x);
y1=zeros(1,xN);   % initiate incoming radiation for ice
y2=zeros(1,xN);   % initiate incoming radiation for cloud
y3=zeros(1,xN);   % initiate outgoing radiation

xsig=5.67E-8;   % Stephen-Boltzmann constant
xRs=340.25;     % insolation constant Q/4

% y1(1:100)=0.1;
% y1(201:251)=0.5;
% y1(100:200)=0.02*((270:0.2:290)-270)+0.1;
% 
% y2(1:50)=0.5;
% y2(225:251)=0.1;
% y2(50:225)=-0.4/35*((260:0.2:295)-260)+0.5;



for xii=1:xN
    y11=-0.000014*abs(x(xii)-288)^3+0.75;    %0.17*sqrt(1-(x(xii)-283)^2/21^2)+0.66;
    if y11<0.6
        y1(xii)=0.6;
    else
        y1(xii)=y11;
    end

    xk=276.645;   % turning point for OLR, x-axis
    yk=0.72471*xRs;    % turning point y-axis
    squ=0.5;
    slope=2.18;

    y31=slope*(x(xii)-xk)+yk;  % will change it later! (-339.647+2.218x)
    
    if x(xii)>xk
        y3(xii)=yk+(4/squ)*slope*(1/(1+exp(squ*(xk-x(xii))))-0.5);  
    else
        y3(xii)=y31;
    end
    %y3(xii)=y3(xii)-0.0000462809917355372*xRs;
    %y3(xii)=y3(xii)-0*0.00001157024793388429*xRs;
    y3(xii)=y3(xii)-0*0.0007*xRs;
    
    % t=3E6;
    % if t<2.5E6
    %     k=0.0007*(t/2.75E6-1);  % CO2 change
    % else
    %     k=0.0007*(t/2.75E6-20/11)/11;
    % end
    % p=0.000092*(0.45*cos(2*pi*t./100000)+0.55*sin(2*pi*t./41000));  % periodic forcing

end


%plot(x,y1,x,y2)
subplot(1,2,1)
plot(x,y1,'LineWidth',1,'Color','#FFA500')

hold on

plot(x,y3/xRs,'LineWidth',1,'Color','#00BFFF')
hold on

greenX=[257.2,291];
greenY=[0.600125,0.749599];
redX=[272.4];
redY=[0.696812];
% greenX=[257.2,281.4,286.2];
% greenY=[0.600125,0.745986,0.74991];
% redX=[272.6,283.4];
% redY=[0.698793,0.748646];
plot(greenX, greenY, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', '#228B22','LineWidth',2);
hold on
plot(redX, redY, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', '#DC143C','LineWidth',2);


hold off
%ylim([0.74,0.755]);

xlim([255, 295]);
%xlabel('T (K)')
%ylabel('radiation (*340.25Wm-2)')
%grid on
legend('R_{in}','R_{out}')
ax=gca;
ax.FontSize = 14;
subplot(1,2,2)
y4=-y1*xRs+y3;
for xii=1:xN
    y5(xii)=sum(y4(1:xii));
end
%y5=y3/xRs-y1;





plot(x,y5,'LineWidth',2,'Color','#191970')
%grid on

ax=gca;
ax.FontSize = 14;
xlim([255,295])


hold on

greenX=[257.2,281.2,286.2];
greenY=[-288.804,340.317,340.303];
redX=[272.6,283.4];
redY=[382.563,340.482];
% greenX=[257.2,291];
% greenY=[-297.617,293.335];
% redX=[272.2];
% redY=[355.688];
plot(greenX, greenY, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', '#228B22','LineWidth',2);
hold on
plot(redX, redY, 'o', 'MarkerSize', 7, 'MarkerEdgeColor', '#DC143C','LineWidth',2);

hold off
%xlabel('T (K)')
%ylabel('potential')



