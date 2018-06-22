
function err = pend_sim(K)
global Mp Mc L Beq Bp kg kt km rm rmp g Ks

Mp  = 0.21;%
Mc  =0.57;
L   =0.3;%
Beq =4.3;
kg  =3.71;
kt  =0.00767;
km  =0.00767;
rm  =2.6;
rmp =6.35*10^-3;


T_sim = 30;
step_time = T_sim/2;

t = 0;
i = 1;% i=0: position ref, i=1: angle reference
k = 1;%k=0 for down 1 for up

if k == 0
    Ks = 0;
    Bp =0.05;
    g = 9.8;
    %% reference trajectory
    T_x = 1;
    a = 0.07;
    b = 2*pi/T_x;
    ofs = 0;
    phase_ofs = 0.1;

    [amp, T_angle, phase] = find_angle(a, T_x, phase_ofs);
    c = amp;
    d = 2*pi/T_angle;
    
    %% linearization trajectory
    Tx_l = 1;
    a_l = 0.05;
    b_l = 2*pi/Tx_l;
    phase_l = 0;
    ofs_l = 0;

    [amp_l, Tangle_l, phase_l] = find_angle(a_l, Tx_l, phase_l);
    c_l = amp_l;
    d_l = 2*pi/Tangle_l;
    %%
    x = ofs + a*sin(b*t+phase*i+phase_ofs);
    x_d = a*b*cos(b*t+phase*i+phase_ofs);

    alpha = pi*k+c*sin(d*t+phase*abs(i-1));
    alpha_d = c*d*cos(d*t+phase*abs(i-1));


    %%
    x_l = a_l*sin(b_l*t+phase_l*i);
    x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

    alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

else
    Ks = 0.77;
    Bp = 0.1;
    g = 9.8;
    %% reference
    T_angle = 0.2;
    c = 0.1;
    d = 2*pi/T_angle;
    [amp, T_x, phase, ofs] = find_position(c, T_angle, k);
    a = amp;
    b = 2*pi/T_x;
    phase_ofs = 0;
%     phase= 0;
    
    %% linearization
    Tangle_l = T_angle;
    c_l = c;
    d_l = 2*pi/Tangle_l;
    [amp_l, Tx_l, phase_l, ofs_l] = find_position(c_l, Tangle_l, k);
    a_l = amp_l;
    b_l = 2*pi/Tx_l;
%     phase_l = 0;
    x = ofs + a*sin(b*t+phase*i+phase_ofs);
    x_d = a*b*cos(b*t+phase*i+phase_ofs);

    alpha = pi*k+c*sin(d*t+phase*abs(i-1));
    alpha_d = c*d*cos(d*t+phase*abs(i-1));


    %%
    x_l = ofs_l + a_l*sin(b_l*t+phase_l*i);
    x_dl = a_l*b_l*cos(b_l*t+phase_l*i);

    alpha_l = pi*k+c_l*sin(d_l*t+phase_l*abs(i-1));
    alpha_dl = c_l*d_l*cos(d_l*t+phase_l*abs(i-1));

    
end

%%
% initial conditions
dalpha = 0.1;%randn*abs(k)/10
dx = 0.02;%randn*abs(k)/10
init = [x_l+dx alpha_l+dalpha x_dl alpha_dl 0 0];
init_l = [dx dalpha 0 0];
ode_step = 0.001;

K_lin = K;
K_org = K;

[t,y] = ode45(@(t,y) linear_system(t,y,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_lin,step_time),0:ode_step:T_sim,init_l);

[t,z] = ode45(@(t,z) original_system(t,z,T_x,a,ofs,phase_ofs,T_angle,i,k,phase...
    ,amp,Tx_l,a_l,Tangle_l,phase_l,amp_l,K_org,step_time),0:ode_step:T_sim,init);


end
