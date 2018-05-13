%% find position from angle reference

function [amp, T_x, phase,ofs] = find_position(c, T_alpha, k)
% clear
% close all
global Mp Mc L Beq Bp kg kt km rm rmp g

% 
% Mp = 0.21;%
% Mc =0.57;
% L =0.3;%
% Beq =4.3;
% Bp =0.03;%
% kg =3.71;
% kt =0.00767;
% km =0.00767;
% rm =2.6;
% rmp =6.35*10^-3;
% g = 9.8;

I = Mp*L^2/3;
% k=1;
% T_alpha = 0.2;
% c = 0.1;

Tsim = 50;

d=2*pi/T_alpha;

fs = 100;
Ls = Tsim*fs;

t = T_alpha/4;

alpha = pi*k+c*sin(d*t);
alpha_d = c*d*cos(d*t);
alpha_dd = -c*(d^2)*sin(d*t);

x_dd = ((I+Mp*L^2)*alpha_dd + Mp*g*L*sin(alpha) + Bp*alpha_d)/(Mp*L*cos(alpha));
x_d = x_dd*T_alpha/(2*pi);

% x = a*sin(b*t)
[t,y] = ode45(@(t,y) findposition(t,y,d,c,k),(0:(Ls-1))/fs,[0 x_d]);
% k=0 for downward k=1 for upward pendulum

x = y(Ls/2+1:end,1);

Ls = length(x);

F = fft(x)/Ls;
F1 = 2*abs(F);
F1(1) = F(1)/2;
F1 = F1(1:Ls/2+1);
[pks,lc] = findpeaks(F1);
[amp, loc] = max(pks);

frq=fs*(0:(Ls/2))/Ls;

% figure
% plot(fs*(0:(Ls/2))/Ls,F1)


freq_estimate = frq(lc(loc));
T_x = 1/freq_estimate;

phase = angle(F(lc(loc))) + pi/2;

% a = amp;
% b = 2*pi/T_x;
% pos = a*sin(b*t+phase);
% alph = pi*k+c*sin(d*t);
% line = y(:,1)-pos;
% grd = (line(end)-line(1))/Tsim;
% 
% figure
% plot(t,y(:,1),'b')
% hold on
% plot(t,alph,'r')
% hold on
% plot(t,pos+grd*t+line(1),'g')
% hold on
% plot(t,line,'m')
% hold on
% plot(t,pos,'r')
% 
% ofs = mean(y(:,1));
ofs = 0;

end
