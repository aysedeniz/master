function dxdt = original_system1(t,z,T_x,cx,phase_x,T_alpha,calpha,phase_alpha,T_xdot,cxdot,phase_xdot,xdot_ofs...
    ,T_alphadot,calphadot,phase_alphadot,T_v,c_v,phase_v,V_ofs,K)
% t
global Mp Mc L Beq Bp kg kt km rm rmp g

% I = 4*Mp*L/3;
I = Mp*L^2/3;

x = z(1);
alpha=z(2);
x_d=z(3);
alpha_d=z(4);

    dx = 2*pi/T_x;
    dxdot = 2*pi/T_xdot;
    x_desired = cx*sin(dx*t+phase_x);
    x_d_desired = xdot_ofs+cxdot*sin(dxdot*t+phase_xdot);

    dalpha = 2*pi/T_alpha;
    dalphadot = 2*pi/T_alphadot;
    
    alpha_desired = pi+calpha*sin(dalpha*t+phase_alpha);
    alpha_d_desired = calphadot*sin(dalphadot*t+phase_alphadot);
    d_v = 2*pi/T_v;
    V_o = V_ofs+c_v*sin(d_v*t+phase_v);



e = -[x-x_desired alpha-alpha_desired x_d-x_d_desired alpha_d-alpha_d_desired]';

V_f = [10 K]*e;
V = V_f +V_o;

denom = (Mc+Mp)*I+Mc*Mp*L^2+Mp^2*L^2*(sin(alpha))^2;

nom_x2 = (Mp^2*L^3+I*Mp*L)*sin(alpha)*alpha_d^2;
nom_x3 = Mp*L*cos(alpha)*Bp*alpha_d + Mp^2*L^2*g*cos(alpha)*sin(alpha);

nom_alpha1 = -(Mc+Mp)*(Mp*L*g*sin(alpha)+Bp*alpha_d);
nom_alpha3 = Mp*L*cos(alpha)*(-Mp*L*sin(alpha)*alpha_d^2+Beq*x_d);

nom_alpha2 = -(-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V/(rm*rmp))*Mp*L*cos(alpha);

nom_x1 = (I+Mp*L^2)*(-Beq*x_d-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V/(rm*rmp));

x_dd = (nom_x1+nom_x2+nom_x3)/denom;
alpha_dd = (nom_alpha1+nom_alpha2+nom_alpha3)/denom;

dxdt = [x_d alpha_d x_dd alpha_dd]';

end