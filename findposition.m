%% find position ref from angle

function dydt = findposition(t,y,d,c,k)

global Mp Mc L Beq Bp kg kt km rm rmp g
I = Mp*L^2/3;

alpha = pi*k+c*sin(d*t);
alpha_d = c*d*cos(d*t);
alpha_dd = -c*(d^2)*sin(d*t);

x_dd = -((I+Mp*L^2)*alpha_dd + Mp*g*L*sin(alpha) + Bp*alpha_d)/(Mp*L*cos(alpha));

dydt = [y(2); x_dd];

end