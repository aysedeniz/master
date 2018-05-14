
function V = openloop_input(a, T_x, amp, T_angle,i,k, phase,phase_ofs, t)

global Mp Mc L Beq Bp kg kt km rm rmp g

% a : amplitude of x
b = 2*pi/T_x;

% x = a*sin(b*t);
x_d = a*b*cos(b*t+phase*i+phase_ofs);
x_dd = -a*b*b*sin(b*t+phase*i+phase_ofs);

c = amp;
d = 2*pi/T_angle;

alpha = pi*k+c*sin(d*t+phase*abs(i-1));
alpha_d = c*d*cos(d*t+phase*abs(i-1));
alpha_dd = -c*d*d*sin(d*t+phase*abs(i-1));

Fc = (Mc+Mp)*x_dd+Beq*x_d + Mp*L*cos(alpha).*alpha_dd-Mp*L*sin(alpha).*(alpha_d.*alpha_d);

V = rm*rmp*(Fc+(kg^2)*kt*km*x_d/(rm*rmp^2))/(kg*kt);

end
