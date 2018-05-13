function dzdt = linear_system(t,y,T_x,a,ofs,phase_ofs,T_angle,i,k,phase,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K,step_time)

x = y(1);
alpha=y(2);
x_d = y(3);
alpha_d = y(4);

    b_l = 2*pi/Tx_l;

    x_l = a_l*sin(b_l*t+phase_l*i);
    x_d_l = a_l*b_l*cos(b_l*t+phase_l*i);

    c_l = amp_l;
    d_l = 2*pi/Tangle_l;
    
    alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_d_l = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));
    
if t>=step_time
    b = 2*pi/T_x;

    x_s = ofs+a*sin(b*t+phase*i+phase_ofs);
    x_d_s = a*b*cos(b*t+phase*i+phase_ofs);

    c = amp;
    d = 2*pi/T_angle;

    alpha_s = pi*k + c*sin(d*t+phase*abs(i-1));
    alpha_d_s = c*d*cos(d*t+phase*abs(i-1));
    
    e = y-[x_s-x_l alpha_s-alpha_l x_d_s-x_d_l alpha_d_s-alpha_d_l]';
    V_o = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase,phase_ofs, t);
    [A,B]=linearized_model(alpha_s,x_d_s,alpha_d_s,V_o);
else
    e = y;
    V_o = openloop_input(a_l, Tx_l, amp_l, Tangle_l,i,k, phase_l,0, t);
    [A,B]=linearized_model(alpha_l,x_d_l,alpha_d_l,V_o);
end
    
%     e = y;
%e = [x-x_desired alpha-alpha_desired x_d-x_d_desired alpha_d-alpha_d_desired]';


V_f = -K*e;

V = V_f;


dzdt = A*y + B*V;


end