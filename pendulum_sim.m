%% Pendulum simulation
close all
clear
clc

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


T_sim = 20;
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
    T_angle = 0.2;
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

K = [0  0  0  0];

[t,z1] = ode45(@(t,z1) original_system(t,z1,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K,step_time),0:ode_step:T_sim,init);

K = state_feedback4(T_angle,c_l,phase_l,T_x,a,i,k)
 
% K = [92.2150  -15.2291    0    0];
% K = [204.2152   58.9257  -27.3738   25.5122];
% -27.3738   25.5122
K = [71.507500 -10.814800 6.076500 2.234900];

K_lin = K;
K_org = K;

[t,y] = ode45(@(t,y) linear_system(t,y,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_lin,step_time),0:ode_step:T_sim,init_l);

[t,z] = ode45(@(t,z) original_system(t,z,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_org,step_time),0:ode_step:T_sim,init);

%%
% xn = z(1:T_sim/(2*ode_step)+1,2);
% alphan = z(T_sim/(2*ode_step)+1:end,2);
% xdot = z(T_sim/(2*ode_step)+1:end,3);
% alphadot = z(T_sim/(2*ode_step)+1:end,4);
% 
% fs = 1/ode_step;
% Ls = length(alphan);
% 
% F = fft(xn)/Ls;
% F1 = 2*abs(F);
% F1(1) = F(1)/2;
% F1 = F1(1:Ls/2+1);
% [pks,lc] = findpeaks(F1);
% [amp, loc] = max(pks);
% % figure
% % plot(fs*(0:(Ls/2))/Ls,F1)
% frq=fs*(0:(Ls/2))/Ls;
% 
% freq_estimate = frq(lc(loc));
% T_x = 1/freq_estimate;
% phase_x = angle(F(lc(loc))) + pi/2;
% 
% cx = amp;
% dx = 2*pi/T_x;
% xn = cx*sin(dx*t+phase_x);
% 
% F = fft(alphan)/Ls;
% F1 = 2*abs(F);
% F1(1) = F(1)/2;
% F1 = F1(1:Ls/2+1);
% [pks,lc] = findpeaks(F1);
% [amp, loc] = max(pks);
% % figure
% % plot(fs*(0:(Ls/2))/Ls,F1)
% frq=fs*(0:(Ls/2))/Ls;
% 
% freq_estimate = frq(lc(loc));
% T_alpha = 1/freq_estimate;
% phase_alpha = angle(F(lc(loc))) + pi/2;
% 
% calpha = amp;
% dalpha = 2*pi/T_alpha;
% alph = pi+calpha*sin(dalpha*t+phase_alpha);
% 
% F = fft(xdot)/Ls;
% F1 = 2*abs(F);
% F1(1) = F(1)/2;
% F1 = F1(1:Ls/2+1);
% [pks,lc] = findpeaks(F1);
% [amp, loc] = max(pks);
% % figure
% % plot(fs*(0:(Ls/2))/Ls,F1)
% frq=fs*(0:(Ls/2))/Ls;
% 
% freq_estimate = frq(lc(loc));
% T_xdot = 1/freq_estimate;
% phase_xdot = angle(F(lc(loc))) + pi/2;
% 
% cxdot = amp;
% dxdot = 2*pi/T_alpha;
% xdot_ofs = mean(xdot);
% xdot = xdot_ofs+cxdot*sin(dxdot*t+phase_xdot);
% 
% F = fft(alphadot)/Ls;
% F1 = 2*abs(F);
% F1(1) = F(1)/2;
% F1 = F1(1:Ls/2+1);
% [pks,lc] = findpeaks(F1);
% [amp, loc] = max(pks);
% % figure
% % plot(fs*(0:(Ls/2))/Ls,F1)
% frq=fs*(0:(Ls/2))/Ls;
% 
% freq_estimate = frq(lc(loc));
% T_alphadot = 1/freq_estimate;
% phase_alphadot = angle(F(lc(loc))) + pi/2;
% 
% calphadot = amp;
% dalphadot = 2*pi/T_alphadot;
% alphadot = calphadot*sin(dalphadot*t+phase_alphadot);

%%

x_l = ofs_l*i+a_l*sin(b_l*t+phase_l*i);
x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

x_step = ofs +a*sin(b*t+phase*i+phase_ofs);
x_dstep = a*b*cos(b*t+phase*i+phase_ofs);

alpha_step = pi*k+c*sin(d*t+phase*abs(i-1));
alpha_dstep = c*d*cos(d*t+phase*abs(i-1));

e(1:step_time/ode_step,1)=-(z(1:step_time/ode_step,1)-x_l(1:step_time/ode_step));
e(1:step_time/ode_step,2)=-(z(1:step_time/ode_step,2)-alpha_l(1:step_time/ode_step));
e(1:step_time/ode_step,3)=-(z(1:step_time/ode_step,3)-x_dl(1:step_time/ode_step));
e(1:step_time/ode_step,4)=-(z(1:step_time/ode_step,4)-alpha_dl(1:step_time/ode_step));

e(1+step_time/ode_step:length(t),1)=-(z(1+step_time/ode_step:length(t),1)-x_step(1+step_time/ode_step:length(t)));
e(1+step_time/ode_step:length(t),2)=-(z(1+step_time/ode_step:length(t),2)-alpha_step(1+step_time/ode_step:length(t)));
e(1+step_time/ode_step:length(t),3)=-(z(1+step_time/ode_step:length(t),3)-x_dstep(1+step_time/ode_step:length(t)));
e(1+step_time/ode_step:length(t),4)=-(z(1+step_time/ode_step:length(t),4)-alpha_dstep(1+step_time/ode_step:length(t)));

V_o(1:step_time/ode_step) = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase_l,0, (0:(step_time/ode_step)-1)*ode_step);
V_o(1+step_time/ode_step:length(t)) = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase,phase_ofs, step_time:ode_step:T_sim);
V_f = e*K';

Vtemp = V_o'+V_f;

V=Vtemp(T_sim/(2*ode_step)+1:end);

V_ofs = mean(V);

fs = 1/ode_step;
Ls = length(V);

F = fft(V)/Ls;
F1 = 2*abs(F);
F1(1) = F(1)/2;
F1 = F1(1:Ls/2+1);
[pks,lc] = findpeaks(F1);
[amp, loc] = max(pks);
% figure
% plot(fs*(0:(Ls/2))/Ls,F1)
frq=fs*(0:(Ls/2))/Ls;

freq_estimate = frq(lc(loc));
T_v = 1/freq_estimate;
phase_v = angle(F(lc(loc))) + pi/2;
c_v = amp;

%% output info
% z_s(:,1) = dx-(z(:,1)-a*sin(b*t));
% z_s(:,2) = dalpha-(z(:,2)-c*sin(d*t+phase));
% 
% intsqr_x = trapz(t,(z(:,1)-a*sin(b*t)).*(z(:,1)-a*sin(b*t)));
% intsqr_a = trapz(t,(z(:,2)-c*sin(d*t+phase)).*(z(:,2)-c*sin(d*t+phase)));
% 
% Sx = stepinfo(z_s(:,1),t);
% Sa = stepinfo(z_s(:,2),t);
% 
% stx = Sx.SettlingTime;
% sta = Sa.SettlingTime;
% 
% osx = abs(Sx.Overshoot);
% osa = abs(Sa.Overshoot);
% 
% fileID = fopen('outputinfo.txt','a');
% 
% fprintf(fileID,'K = [%f %f %f %f], initial cond = [%f %f %f %f] stx = %f, sta = %f, osx = %f, osa = %f, intsqrx = %f, intsqra = %f \r\n',...
% K,init,stx,sta,osx,osa,intsqr_x,intsqr_a);
% 
% fclose(fileID);
%%

x(1:step_time/ode_step,1) = x_l(1:step_time/ode_step);
x(step_time/ode_step+1:length(t)) = x_step(step_time/ode_step+1:length(t));
xy(1:step_time/ode_step,1) = x_l(1:step_time/ode_step);
xy(step_time/ode_step+1:length(t)) = x_step(step_time/ode_step+1:length(t))-ofs+ ofs_l*i;


alpha(1:step_time/ode_step,1) = alpha_l(1:step_time/ode_step);
alpha(step_time/ode_step+1:length(t)) = alpha_step(step_time/ode_step+1:length(t));

figure(10)
plot(t,V_o,'b');
hold on
plot(t,V_f,'m');
hold on
plot(t,V_o'+V_f,'g')
legend('openloop input','feedback','inp')
grid on

figure(20)
subplot(2,1,1)
plot(t,z(:,1))
hold on
plot(t,x,'r')
title('position')
hold on
plot(t,y(:,1)+xy,'m')
% hold on
% plot(t,V_f/10,'c')
legend('original','reference','linear')
grid on

subplot(2,1,2)
plot(t,z(:,2))
hold on
plot(t,alpha,'r')
hold on
plot(t,y(:,2)+alpha,'m')
title('angle')
grid on

x(step_time/ode_step+1:length(t)) = x_step(step_time/ode_step+1:length(t))-ofs*abs(k-1);

% y(step_time/ode_step+1:length(t),1) = y(step_time/ode_step+1:length(t),1)+ofs;
% v(step_time/ode_step+1:length(t),1) = v(step_time/ode_step+1:length(t),1)+ofs;

% ref(1:step_time/ode_step,1) = x(1:step_time*100,1)-x(1:step_time/ode_step,1);
% ref(step_time/ode_step+1:length(t),1) = ofs*ones(step_time/ode_step+1:length(t));

figure(30)
subplot(2,2,1)
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

subplot(2,2,2)
plot(t,z(:,2)-alpha,'b')
xlabel('t(sn)')
ylabel('angle(rad)')
grid on
hold on
plot(t,y(:,2),'m')
grid on
% hold on
% plot(t,-V_f/100,'c')

subplot(2,2,3)
plot(t,z(:,3)-x_dl,'b')
ylabel('velocity')
grid on
% hold on
% plot(t,-V_f/100,'c')
% hold on
title('error')
legend('original','linear')
grid on

subplot(2,2,4)
plot(t,z(:,4)-alpha_dl,'b')
ylabel('angular velocity')
grid on
% hold on
% plot(t,-V_f/100,'c')
% hold on
title('error')
legend('original','linear')
grid on

s=1;
figure(31)
subplot(2,1,1)
plot(t,z(:,1),'b','LineWidth',s)
set(gca,'FontSize',12)
ylabel('pozisyon(m)','FontSize',14)
hold on
plot(t,z1(:,1),'m','LineWidth',s)
grid on
legend('geri beslemeli','geri beslemesiz')


subplot(2,1,2)
plot(t,z(:,2),'b','LineWidth',s)
xlabel('t(sn)','FontSize',14)
ylabel('açý(rad)','FontSize',14)
grid on
hold on
plot(t,z1(:,2),'m','LineWidth',s)
grid on

figure(40)
subplot(2,1,1)
plot(t,z(:,3))
title('translational velocity')
grid on
% hold on
% plot(t,y(:,3),'m')
% hold on
% plot(t,xdot,'g')

subplot(2,1,2)
plot(t,z(:,4))
title('angular velocity')
grid on
% hold on
% plot(t,y(:,4),'m')
% hold on
% plot(t,alphadot,'g')

figure(60)
if k == 0
    plot(z(1:step_time/ode_step,1),z(1:step_time/ode_step,2),'m')
    hold on
    plot(z(step_time/ode_step+1:end,1),z(step_time/ode_step+1:end,2),'g')
else
    plot(z(:,1),z(:,2),'m')
end

hold on
plot(x(1:step_time/ode_step),alpha(1:step_time/ode_step),'b','LineWidth',1.2)
hold on
plot(x(step_time/ode_step+1:length(t)),alpha(step_time/ode_step+1:length(t)),'b','LineWidth',1.2)
hold on
% plot(z(1,1),z(2,2),'>m','LineWidth',2)
xlabel('pozisyon(m)','FontSize',14)
ylabel('açý(rad)','FontSize',14)
% legend('sistem çýktýsý','referans yörüngesi')
P=1000;
quiver(z(1,1),z(1,2),z(1,3)/100,z(1,4)/100,'m','MaxHeadSize',3000,'LineWidth',1.2)
hold on
quiver(z(350,1),z(350,2),z(350,3)/100,z(350,4)/100,'m','MaxHeadSize',2000,'LineWidth',1.2)
hold on
quiver(z(P,1),z(P,2),z(P,3)/100,z(P,4)/100,'m','MaxHeadSize',2000,'LineWidth',1.2)
grid on

figure
subplot(2,1,1)
plot(t,z(:,5))
subplot(2,1,2)
plot(t,z(:,6))
