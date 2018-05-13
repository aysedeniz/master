function dvdt = linear_system_approx(t,v,T_x,a,T_angle,i,k,phase,amp,K,A,B)

% x = y(1);
% alpha = y(2);
% x_d = y(3);
% alpha_d = y(4);

b = 2*pi/T_x;
% 
% x_desired = a*sin(b*t+phase*i);
% x_d_desired = a*b*cos(b*t+phase*i);

c = amp;
d = 2*pi/T_angle;

% alpha_desired = pi*k+c*sin(d*t+phase*abs(i-1));
% alpha_d_desired = c*d*cos(d*t+phase*abs(i-1));

%e = [x-x_desired alpha-alpha_desired x_d-x_d_desired alpha_d-alpha_d_desired]';
e = v;

U = -K*e;
Atot = zeros(4,4);
Btot = zeros(4,1);
for j=-2:2
    Atot = Atot + A(:,:,j+3)*exp(1i*b*t*j);
    Btot = Btot + B(:,:,j+3)*exp(1i*b*t*j);
end

dvdt = Atot*v + Btot*U;


end