%%
function [amp, T_angle, phase] = find_angle(a, T_x, p)

global Mp Mc L Beq Bp kg kt km rm rmp g
I = Mp*L^2/3;
% 
% T_x = 0.2;
% a = 0.05;
% p = 0.5;

b=2*pi/T_x;

fs = 100;
Ls = 20*fs;
t= 0;

[t,y] = ode45(@(t,y) findangle(t,y,b,a,p),(0:(Ls-1))/fs,[0 0]);

alpha = y(Ls/2+1:end,1);

Ls = length(alpha);

F = fft(alpha)/Ls;
F1 = 2*abs(F);
F1(1) = F(1)/2;
F1 = F1(1:Ls/2+1);
[pks,lc] = findpeaks(F1);
[amp, loc] = max(pks);

frq=fs*(0:(Ls/2))/Ls;

% F1(lc(loc))=0;

% figure
% plot(fs*(0:(Ls/2))/Ls,F1)


freq_estimate = frq(lc(loc));
T_angle = 1/freq_estimate

phase = angle(F(lc(loc))) + pi/2

c = amp
d = 2*pi/T_angle;
x = a*sin(b*t+p);
alph = c*sin(d*t+phase);

figure
plot(t,y(:,1),'b','LineWidth',1.5)
hold on
plot(t,alph,'m','LineWidth',1.5)
legend('numerical integration output','fitted sinusoid')
grid on


end

