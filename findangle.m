function dydt = findangle(t,y,b,a,p)

global Mp Mc L Beq Bp kg kt km rm rmp g
I = Mp*L^2/3;

alpha = y(1);
alpha_d = y(2);
x = -a*(b^2)*sin(b*t+p);
alpha_dd = (-Bp*alpha_d-Mp*g*L*sin(alpha)-Mp*L*cos(alpha)*x)/(I+Mp*L^2);
dydt = [y(2); alpha_dd];

end