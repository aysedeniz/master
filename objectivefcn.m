function f = objectivefcn(K)
global Mp Mc L Beq Bp kg kt km rm rmp g

Mp = 0.21;%
Mc =0.57;
L =0.3;%
Beq =4.3;
Bp =0.03;%
kg =3.71;
kt =0.00767;
km =0.00767;
rm =2.6;
rmp =6.35*10^-3;
g = 9.8;

T_sim = 30;

t = 0;
i = 0;% i=0: position ref, i=1: angle reference
k = 0;%k=0 for down 1 for up
%%
T_x = 0.2;
a = 0.05;
b = 2*pi/T_x;

[amp, T_angle, phase] = find_angle(a, T_x);
c = amp;
d = 2*pi/T_angle;

%%
% T_angle = 2;
% c = 0.2;
% d = 2*pi/T_angle;
% [amp, T_x, phase] = find_position(c, T_angle, k);
% a = amp;
% b = 2*pi/T_angle;
%%
x = a*sin(b*t+phase*i);
x_d = a*b*cos(b*t+phase*i);

alpha = pi*k+c*sin(d*t+phase*abs(i-1));
alpha_d = c*d*cos(d*t+phase*abs(i-1));
dx = 0.1;
dalpha = 0.1;

% [t,z] = ode45(@(t,z) original_system(t,z,T_x,a,T_angle,i,k,phase,amp,K),0:0.01:T_sim,[x+dx alpha+dalpha x_d alpha_d]);
% 
% 
% z(:,1) = dx-(z(:,1)-a*sin(b*t));
% z(:,2) = dalpha-(z(:,2)-c*sin(d*t+phase));

[t,z] = ode45(@(t,z) linear_system(t,z,T_x,a,T_angle,i,k,phase,amp,K),0:0.01:T_sim,[dx dalpha 0 0]);

z(:,1) = dx-z(:,1);
z(:,2) = dalpha-z(:,2);

Sx = stepinfo(z(:,1),t);
Sa = stepinfo(z(:,2),t);

stx = Sx.SettlingTime
sta = Sa.SettlingTime

osx = abs(Sx.Overshoot)
osa = abs(Sa.Overshoot)

f = (200*(stx + sta) + osx + osa)/1000;
K

end
