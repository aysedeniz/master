% Linearization
function [A,B]=linearized_model_3state(alpha,x_d,alpha_d,V_x)
global Mp Mc L Beq Bp kg kt km rm rmp g
I = Mp*L^2/3;

nom_x1 = (I+Mp*L^2)*(-Beq*x_d-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V_x/(rm*rmp));
nom_x2 = (Mp^2*L^3+I*Mp*L)*sin(alpha)*alpha_d^2;
nom_x3 = Mp*L*cos(alpha)*Bp*alpha_d + Mp^2*L^2*g*cos(alpha)*sin(alpha);

nom_alpha1 = -(Mc+Mp)*(Mp*L*g*sin(alpha)+Bp*alpha_d);
nom_alpha2 = -(-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V_x/(rm*rmp))*Mp*L*cos(alpha);
nom_alpha3 = Mp*L*cos(alpha)*(-Mp*L*sin(alpha)*alpha_d^2+Beq*x_d);

denom = (Mc+Mp)*I+Mc*Mp*L^2+Mp^2*L^2*(sin(alpha))^2;


xdd1 = (((Mp^2*L^3+I*Mp*L)*cos(alpha)*alpha_d^2 - Mp*L*sin(alpha)*Bp*alpha_d + Mp^2*L^2*g*cos(2*alpha))/denom) - Mp^2*L^2*sin(2*alpha)*(nom_x1+nom_x2+nom_x3);
xdd2 = (I+Mp*L^2)*(-Beq-kg^2*kt*km/(rm*rmp^2))/denom;
xdd3 = (2*(Mp^2*L^3+I*Mp*L)*sin(alpha)*alpha_d + Mp*L*cos(alpha)*Bp)/denom;

L1 = -Mp*L*sin(alpha)*(-Mp*L*sin(alpha)*alpha_d^2+Beq*x_d);
L2 = Mp*L*cos(alpha)*(-Mp*L*cos(alpha)*alpha_d^2);

alphadd1 = ((-(Mc+Mp)*(Mp*L*g*cos(alpha))+(-kg^2*kt*km*x_d/(rm*rmp^2)+kg*kt*V_x/(rm*rmp))*Mp*L*sin(alpha)+L1+L2)/denom)- Mp^2*L^2*sin(2*alpha)*(nom_alpha1+nom_alpha2+nom_alpha3);
alphadd2 = ((kg^2*kt*km/(rm*rmp^2))*Mp*L*cos(alpha) + Mp*L*cos(alpha)*Beq*x_d)/denom;
alphadd3 = (-(Mc+Mp)*Bp + 2*Mp*L*cos(alpha)*(-Mp*L*sin(alpha)*alpha_d))/denom;

B1 = 0;
B2 = ((I+Mp*L^2)*kg*kt/(rm*rmp))/denom;
B3 = -(kg*kt/(rm*rmp))*Mp*L*cos(alpha)/denom;


A = [0 0 1; xdd1 xdd2 xdd3; alphadd1 alphadd2 alphadd3];
B = [B1;B2;B3];


