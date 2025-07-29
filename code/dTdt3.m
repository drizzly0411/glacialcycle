function [F,T_high,T_unstable,T_low]=dTdt3(t,T)
    C=4.2E8;
    xRs=340.25;
    y11=@(x) -0.000014*abs(x-288)^3+0.75;
    if y11(T)<0.6
        y1=0.6;
    else
        y1=y11(T);
    end
    xk=276.645;   % turning point for OLR, x-axis
    yk=0.72471*xRs;    % turning point y-axis
    squ=0.5;
    slope=2.18;

    y31=@(x) slope*(x-xk)+yk;   % see Koll et al., moved
    y32=@(x) yk+(4/squ)*slope*(1/(1+exp(squ*(xk-x)))-0.5);
    


    if T>xk
        y3=y32(T);
    else
        y3=y31(T);
    end
    
    if t<2.5E6
        k=0.0007*(t/2.75E6-1);  % CO2 change
    else
        k=0.0007*(t/2.75E6-20/11)/11;
    end
    
            q=0*0.00002*(283-T);     % CO2 change feedback (closed)
    p=0.000092*(0.45*cos(2*pi*t./100000)+0.55*sin(2*pi*t./41000));  % periodic forcing

    F=((1+p)*y1*xRs-(y3+(k+q)*xRs))/C;     % dT/dt=F, at particular T and t
    %F=-0.14*(T-291.057)/C;
    %F=-0.0274421*(T-286.302)/C;
    
    %-----------------------------------------
    % % find the zero points
    % FF=@(x) ((1+p)*y11(x)*340.25-(y32(x)+(k+q)*340.25))*10000;
    % %rootgrid=
    % T_unstable=fzero(FF,281);
    % T_high=fzero(FF,283);
    % %FF2=@(x) ((1+p)*y11(x)*340.25-(y31(x)+k*340.25))*10000;
    % T_low=fzero(FF,280);
    T_high=0;
    T_unstable=0;
    T_low=0;    % temporary value
    

end