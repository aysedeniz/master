
function dxdt = original_system(t,z,T_x,a,ofs,phase_ofs,T_angle,i,k , phase, amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K,step_time)
% t
global Mp Mc L Beq Bp kg kt km rm rmp g Ks

% I = 4*Mp*L/3;
I = Mp*L^2/3;

x=z(1);
alpha=z(2);
x_d=z(3);
alpha_d=z(4);
if t>=step_time

    b = 2*pi/T_x;

    x_desired = ofs+a*sin(b*t+phase*i+phase_ofs);
    x_d_desired = a*b*cos(b*t+phase*i+phase_ofs);

    c = amp;
    d = 2*pi/T_angle;

    alpha_desired = pi*k + c*sin(d*t+phase*abs(i-1));
    alpha_d_desired = c*d*cos(d*t+phase*abs(i-1));
    
%     c_l = amp_l;
%     d_l = 2*pi/Tangle_l;
    
%     alpha_desired = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
%     alpha_d_desired = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));
    V_o = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase,phase_ofs, t);
else
    b_l = 2*pi/Tx_l;

    x_desired = a_l*sin(b_l*t+phase_l*i);
    x_d_desired = a_l*b_l*cos(b_l*t+phase_l*i);

    c_l = amp_l;
    d_l = 2*pi/Tangle_l;
    
    alpha_desired = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_d_desired = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));
    V_o = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase_l,0, t);

end

e = -[x-x_desired alpha-alpha_desired x_d-x_d_desired alpha_d-alpha_d_desired]';

V_f = K*e + z(5) + z(6);
V = V_o + V_f;

denom = (Mc+Mp)*I+Mc*Mp*L^2+Mp^2*L^2*(sin(alpha))^2;

nom_x2 = (Mp^2*L^3+I*Mp*L)*sin(alpha)*alpha_d^2;
nom_x3 = Mp*L*cos(alpha)*(Bp*alpha_d+Ks*(alpha-pi)) + Mp^2*L^2*g*cos(alpha)*sin(alpha);

nom_alpha1 = -(Mc+Mp)*(Mp*L*g*sin(alpha)+Bp*alpha_d+Ks*(alpha-pi));
nom_alpha3 = Mp*L*cos(alpha)*(-Mp*L*sin(alpha)*alpha_d^2+Beq*x_d);

nom_alpha2 = -(-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V/(rm*rmp))*Mp*L*cos(alpha);

nom_x1 = (I+Mp*L^2)*(-Beq*x_d-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V/(rm*rmp));

x_dd = (nom_x1+nom_x2+nom_x3)/denom;
alpha_dd = (nom_alpha1+nom_alpha2+nom_alpha3)/denom;
Kix = 0;
Kia = 0;
eix = Kix*e(1);
eia = Kia*e(2);
dxdt = [x_d alpha_d x_dd alpha_dd eix eia]';

end