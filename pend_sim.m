
function [err, conv_err] = pend_sim(K,T)
global Mp Mc L Beq Bp kg kt km rm rmp g Ks

Mp  = 0.21;%
Mc  =0.57;
L   =0.3;%
Beq =4.3;
kg  =3.71;
kt  =0.00767;
km  =0.00767;
rm  =2.6;
rmp =6.35*10^-3;

if nargin<2
    T = 0.2;
end
T_sim = 30;
step_time = T_sim/2;

t = 0;
i = 1;% i=0: position ref, i=1: angle reference
k = 1;%k=0 for down 1 for up

if k == 0
    Ks = 0;
    Bp =0.05;
    g = 9.8;
    %% reference trajectory
    T_x = 1;
    a = 0.07;
    b = 2*pi/T_x;
    ofs = 0;
    phase_ofs = 0.1;

    [amp, T_angle, phase] = find_angle(a, T_x, phase_ofs);
    c = amp;
    d = 2*pi/T_angle;
    
    %% linearization trajectory
    Tx_l = 1;
    a_l = 0.05;
    b_l = 2*pi/Tx_l;
    phase_l = 0;
    ofs_l = 0;

    [amp_l, Tangle_l, phase_l] = find_angle(a_l, Tx_l, phase_l);
    c_l = amp_l;
    d_l = 2*pi/Tangle_l;
    %%
    x = ofs + a*sin(b*t+phase*i+phase_ofs);
    x_d = a*b*cos(b*t+phase*i+phase_ofs);

    alpha = pi*k+c*sin(d*t+phase*abs(i-1));
    alpha_d = c*d*cos(d*t+phase*abs(i-1));


    %%
    x_l = a_l*sin(b_l*t+phase_l*i);
    x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

    alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

else
    Ks = 0.77;
    Bp = 0.1;
    g = 9.8;
    %% reference
    T_angle = T;
    c = 0.1;
    d = 2*pi/T_angle;
    [amp, T_x, phase, ofs] = find_position(c, T_angle, k);
    a = amp;
    b = 2*pi/T_x;
    phase_ofs = 0;
%     phase= 0;
    
    %% linearization
    Tangle_l = T_angle;
    c_l = c;
    d_l = 2*pi/Tangle_l;
    [amp_l, Tx_l, phase_l, ofs_l] = find_position(c_l, Tangle_l, k);
    a_l = amp_l;
    b_l = 2*pi/Tx_l;
%     phase_l = 0;
    x = ofs + a*sin(b*t+phase*i+phase_ofs);
    x_d = a*b*cos(b*t+phase*i+phase_ofs);

    alpha = pi*k+c*sin(d*t+phase*abs(i-1));
    alpha_d = c*d*cos(d*t+phase*abs(i-1));


    %%
    x_l = ofs_l + a_l*sin(b_l*t+phase_l*i);
    x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

    alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

    
end

%%
% initial conditions
dalpha = 0.1;%randn*abs(k)/10
dx = 0.02;%randn*abs(k)/10
init = [x_l+dx alpha_l+dalpha x_dl alpha_dl 0 0];
init_l = [dx dalpha 0 0];
ode_step = 0.001;

K_lin = K
K_org = K;

opt    = odeset('Events', @event_unstable);

[t,z,~,~,ie] = ode45(@(t,z) original_system(t,z,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_org,step_time),0:ode_step:T_sim,init,opt);
ie
if isempty(ie)
    [t,y] = ode45(@(t,y) linear_system(t,y,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
        ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_lin,step_time),0:ode_step:T_sim,init_l);

    x = z(:,1)-(ofs + a*sin(b*t+phase*i+phase_ofs));
    for j = 1:(length(z(:,1))-200)
        x_smooth(j) = mean(x(j:j+200));
    end

    diff_x_smooth = x_smooth'-y(1:length(x_smooth),1);
    diff_x = x(1:length(x_smooth))-y(1:length(x_smooth),1);
    alpha = pi*k+c*sin(d*t+phase*abs(i-1));
    for j = 1:(length(z(:,2))-200)
        alpha_smooth(j) = mean(alpha(j:j+200));
    end

    diff_alpha_smooth = alpha_smooth'-y(1:length(alpha_smooth),2);
    diff_alpha = alpha(1:length(alpha_smooth))-y(1:length(alpha_smooth),2);
    normal = sum(abs(y(1:length(alpha_smooth),1))+abs(y(1:length(alpha_smooth),2)));
    conv_err = (sum(abs(diff_x_smooth))+sum(abs(diff_alpha_smooth)))/normal
    err = (sum(abs(diff_x)+sum(abs(diff_alpha))))/normal
else
    conv_err = 10000
    err = 10000
end
x_l = ofs_l*i+a_l*sin(b_l*t+phase_l*i);
x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

x_step = ofs +a*sin(b*t+phase*i+phase_ofs);
x_dstep = a*b*cos(b*t+phase*i+phase_ofs);

alpha_step = pi*k+c*sin(d*t+phase*abs(i-1));
alpha_dstep = c*d*cos(d*t+phase*abs(i-1));

x(1:step_time/ode_step,1) = x_l(1:step_time/ode_step);
x(step_time/ode_step+1:length(t)) = x_step(step_time/ode_step+1:length(t));
xy(1:step_time/ode_step,1) = x_l(1:step_time/ode_step);
xy(step_time/ode_step+1:length(t)) = x_step(step_time/ode_step+1:length(t))-ofs+ ofs_l*i;


alpha(1:step_time/ode_step,1) = alpha_l(1:step_time/ode_step);
alpha(step_time/ode_step+1:length(t)) = alpha_step(step_time/ode_step+1:length(t));

if nargin == 2
%     figure
%     subplot(2,1,1)
%     plot(t,z(:,1))
%     hold on
%     plot(t,x,'r')
%     title('position')
%     hold on
%     plot(t,y(:,1)+xy,'m')
%     % hold on
%     % plot(t,V_f/10,'c')
%     legend('original','reference','linear')
%     grid on
% 
%     subplot(2,1,2)
%     plot(t,z(:,2))
%     hold on
%     plot(t,alpha,'r')
%     hold on
%     plot(t,y(:,2)+alpha,'m')
%     title('angle')
%     grid on
figure
subplot(2,1,1)
plot(t,z(:,1)-x,'b')
ylabel('position(m)')
hold on
plot(t,y(:,1)-(x-x_l),'m')
grid on
% hold on
% plot(t,-V_f/100,'c')
% hold on
title('error')
legend('original','linear')
grid on

subplot(2,1,2)
plot(t,z(:,2)-alpha,'b')
xlabel('t(sn)')
ylabel('angle(rad)')
grid on
hold on
plot(t,y(:,2),'m')
grid on
end
end
