clear; clc; close all;
set(0,'DefaultFigureWindowStyle','normal') %'normal' 'docked'
set(0,'defaulttextInterpreter','latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultfigurecolor',[1 1 1])

[A, B, Rxx, Ruu, Pf, tf]=getParams([]);

interval=tf:-0.01:0.0; %Note that we need to integrate backwards. See also https://groups.google.com/g/comp.soft-sys.matlab/c/17Q5OpbueDQ/m/sV27KeBeqQgJ

%%%%%%%%%%%%%%%%%%%% SOLVE THE DIFFERENTIAL RICATTI EQUATION
Pf_flattened=[Pf(1,1) ; Pf(1,2) ; Pf(2,2)];
[t,P_flattened] = ode45(@odefun,interval,Pf_flattened);
K=[];
for i=1:size(P_flattened,1)
   P_t=[P_flattened(i,1), P_flattened(i,2); P_flattened(i,2) P_flattened(i,3)];
   K=[K;
      inv(Ruu)*B'*P_t ];
end

%%%%%%%%%%%%%%%%%%%% SOLVE THE CONTROL ALGEBRAIC RICATTI EQUATION
Q = Rxx;
R = Ruu;
S = []; %Default is 0
E = []; %Default is I
G = []; %Default is 0

[P_ss,K_ss] = icare(A,B,Q,R,S,E,G);

A_closed_loop_ss=(A-B*inv(Ruu)*B'*P_ss);


%%%%%%%%%%%%%%%%%%%%
%PLOTS
figure; hold on;
plot(t,P_flattened,'o')
o=ones(size(t));
plot(t,P_ss(1,1)./o,'LineWidth',2.0)
plot(t,P_ss(1,2)./o,'LineWidth',2.0)
plot(t,P_ss(2,2)./o,'LineWidth',2.0)

xlabel('t (s)')
legend('P(0,0)','P(0,1)','P(1,1)','Pss(0,0)','Pss(0,1)','Pss(1,1)','Location','northeastoutside')
title('Matrix P =[P(0,0) P(0,1) ; P(0,1) P(1,1)]')

%PLOTS
figure; hold on;
plot(t,K,'o')
o=ones(size(t));
plot(t,K_ss(1)./o,'LineWidth',2.0)
plot(t,K_ss(2)./o,'LineWidth',2.0)

xlabel('t (s)')
legend('K(0,0)','K(0,1)','Kss(0,0)','Kss(0,1)','Location','northeastoutside')
title('Matrix K =[K(0,0) K(0,1]')



function result = odefun(t,y)

    [A, B, Rxx, Ruu, Pf, tf]=getParams(t);

    P=[y(1), y(2); y(2) y(3)];

    P_dot=-(P*A +A'*P+Rxx-P*B*inv(Ruu)*B'*P);

    result=[P_dot(1,1) ; P_dot(1,2) ; P_dot(2,2)];
end

function [A, B, Rxx, Ruu, Pf, tf]=getParams(t);
    
    tf=10;

    A=[0 1;
       0 0];
    
    B=[0;
       1];
    
    Rxx=[2 0;
         0 1];
    
    Ruu=[1];

    Pf=[0.7 0.5;
    0.5 1.0];
    
end

