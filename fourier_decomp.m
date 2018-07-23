%% fourier decomposition
clear
close all
clc

global Mp Mc L Beq Bp kg kt km rm rmp g Ks

Mp = 0.21;
Mc =0.57;
L =0.3;
Beq =4.3;
Bp =0.1;
kg =3.71;
kt =0.00767;
km =0.00767;
rm =2.6;
rmp =6.35*10^-3;
g = 9.8;
Ks = 0.77;

I = Mp*L^2/3;

i = 1;
k = 1;

Fs = 400;
Ls = Fs*10;

T_angle = 0.2;
c = 0.1;
d = 2*pi/T_angle;
[amp, T_x, phase] = find_position(c, T_angle, k);
a = amp;
b = 2*pi/T_x;
ofs = 0;
phase_ofs = 0;
T_x =0.2;

% T_x = 1;
% a = 0.05;
% b = 2*pi/T_x;
% p=0;
% 
% [amp, T_angle, phase] = find_angle(a, T_x, p);
% c = amp;
% d = 2*pi/T_angle;

% t = (0:(Ls-1))/Fs;
% x_d = a*b*cos(b*t+phase*i);
% 
% alpha = pi*k+c*sin(d*t+phase*abs(i-1));
% alpha_d = c*d*cos(d*t+phase*abs(i-1));
%    figure
%    plot(t,alpha);
%    figure
%    plot(t,alpha_d);
%    figure
%    plot(t,x_d);
%    
p = 0;
for t = (0:(Ls-1))/Fs  
   p = p+1;
   alpha(p) = pi*k+c*sin(d*t+phase*abs(i-1));
   alpha_d(p) = c*d*cos(d*t+phase*abs(i-1));
   x_d(p) = a*b*cos(b*t+phase*i);
   V(p) = openloop_input(a, T_x, c, T_angle,i,k, phase,0, t);
   [A,B]=linearized_model(alpha(p),x_d(p),alpha_d(p),V(p));
    for k=0:1
        for j=1:3
            A_n(j+(k*3),p) = A(3+k,j+1);
            if j==1
                B_n(k+1,p) = B(k+3);
            end
        end
    end
end

p=(0:Ls-1)/Fs;
figure
plot(p,alpha)
figure
plot(p,alpha_d)
figure
plot(p,x_d)
figure
plot(p,V)
figure
plot(p,A_n(1,:))


%     figure(1)
% for k=1:6
%     subplot(2,3,k)
%     plot((0:Ls-1)/Fs,A_n(k,:))
% end
% 
%     figure(2)
% for k=1:2
%     subplot(2,1,k)
%     plot((0:Ls-1)/Fs,B_n(k,:))
% end

%     figure(3)
for k=1:6
    f(k,:) = fft(A_n(k,:))/Ls;
    f1(k,:) = 2*abs(f(k,:));
    f1(k,1) = f1(k,1)/2;
    f2(k,:) = f1(k,1:Ls/2+1);
    
%     subplot(2,3,k)
%     plot(Fs*(0:(Ls/2))/Ls,f2(k,:))
end

%     figure(4)
for k=1:2
    fb(k,:) = fft(B_n(k,:))/Ls;
    fb1(k,:) = 2*abs(fb(k,:));
    fb1(k,1) = fb1(k,1)/2;
    fb2(k,:) = fb1(k,1:Ls/2+1);
    
%     subplot(2,1,k)
%     plot(Fs*(0:(Ls/2))/Ls,fb2(k,:))
end

loc(1) = (Ls/(T_x*Fs))+1;
loc(2) = 2*(Ls/(T_x*Fs))+1;
loc(3) = 3*(Ls/(T_x*Fs))+1;
loc(4) = 4*(Ls/(T_x*Fs))+1;
loc(5) = 5*(Ls/(T_x*Fs))+1;



for k=1:6
    a1(k) = f2(k,loc(1))*exp(angle(f(k,loc(1)))*1i)/2;
    a_1(k) = f2(k,loc(1))*exp(-angle(f(k,loc(1)))*1i)/2;
end


for k=1:6
    a2(k) = f2(k,loc(2))*exp(angle(f(k,loc(2)))*1i)/2;
    a_2(k) = f2(k,loc(2))*exp(-angle(f(k,loc(2)))*1i)/2;
end

for k=1:6
    a3(k) = f2(k,loc(3))*exp(angle(f(k,loc(3)))*1i)/2;
    a_3(k) = f2(k,loc(3))*exp(-angle(f(k,loc(3)))*1i)/2;
end

for k=1:6
    a4(k) = f2(k,loc(4))*exp(angle(f(k,loc(4)))*1i)/2;
    a_4(k) = f2(k,loc(4))*exp(-angle(f(k,loc(4)))*1i)/2;
end

for k=1:6
    a5(k) = f2(k,loc(5))*exp(angle(f(k,loc(5)))*1i)/2;
    a_5(k) = f2(k,loc(5))*exp(-angle(f(k,loc(5)))*1i)/2;
end

A0 = [0     0       1      0;
      0     0       0      1;
      0 f(1,1) f(2,1) f(3,1);
      0 f(4,1) f(5,1) f(6,1)];
  
B0 = [0; 0; fb(1,1); fb(2,1)];

A1 = [0     0       0      0;
      0     0       0      0;
      0     a1(1)   a1(2)  a1(3);
      0     a1(4)   a1(5)  a1(6)];
  
A_1 = [0     0       0      0;
      0     0       0      0;
      0     a_1(1)   a_1(2)  a_1(3);
      0     a_1(4)   a_1(5)  a_1(6)];
  
B1 = [0; 0; fb2(1,loc(1))*exp(angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(angle(fb(2,loc(1)))*1i)/2];
B_1 = [0; 0; fb2(1,loc(1))*exp(-angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(-angle(fb(2,loc(1)))*1i)/2];


A2 = [0     0       0      0;
      0     0       0      0;
      0     a2(1)   a2(2)  a2(3);
      0     a2(4)   a2(5)  a2(6)];
  
A_2 = [0     0       0      0;
      0     0       0      0;
      0     a_2(1)   a_2(2)  a_2(3);
      0     a_2(4)   a_2(5)  a_2(6)];

B2 = [0; 0; fb2(1,loc(2))*exp(angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(angle(fb(2,loc(2)))*1i)/2];
B_2 = [0; 0; fb2(1,loc(2))*exp(-angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(-angle(fb(2,loc(2)))*1i)/2];

A3 = [0     0       0      0;
      0     0       0      0;
      0     a3(1)   a3(2)  a3(3);
      0     a3(4)   a3(5)  a3(6)];
  
A_3 = [0     0       0      0;
      0     0       0      0;
      0     a_3(1)   a_3(2)  a_3(3);
      0     a_3(4)   a_3(5)  a_3(6)];

B3 = [0; 0; fb2(1,loc(3))*exp(angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(angle(fb(2,loc(3)))*1i)/2];
B_3 = [0; 0; fb2(1,loc(3))*exp(-angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(-angle(fb(2,loc(3)))*1i)/2];

A4 = [0     0       0      0;
      0     0       0      0;
      0     a4(1)   a4(2)  a4(3);
      0     a4(4)   a4(5)  a4(6)];
  
A_4 = [0     0       0      0;
      0     0       0      0;
      0     a_4(1)   a_4(2)  a_4(3);
      0     a_4(4)   a_4(5)  a_4(6)];

B4 = [0; 0; fb2(1,loc(4))*exp(angle(fb(1,loc(4)))*1i)/2;...
    fb2(2,loc(4))*exp(angle(fb(2,loc(4)))*1i)/2];
B_4 = [0; 0; fb2(1,loc(4))*exp(-angle(fb(1,loc(4)))*1i)/2;...
    fb2(2,loc(4))*exp(-angle(fb(2,loc(4)))*1i)/2];

A5 = [0     0       0      0;
      0     0       0      0;
      0     a5(1)   a5(2)  a5(3);
      0     a5(4)   a5(5)  a5(6)];
  
A_5 = [0     0       0      0;
      0     0       0      0;
      0     a_5(1)   a_5(2)  a_5(3);
      0     a_5(4)   a_5(5)  a_5(6)];

B5 = [0; 0; fb2(1,loc(5))*exp(angle(fb(1,loc(5)))*1i)/2;...
    fb2(2,loc(5))*exp(angle(fb(2,loc(5)))*1i)/2];
B_5 = [0; 0; fb2(1,loc(5))*exp(-angle(fb(1,loc(5)))*1i)/2;...
    fb2(2,loc(5))*exp(-angle(fb(2,loc(5)))*1i)/2];


    
 A_t =   [A0 A_1 A_2 A_3 A_4;
          A1 A0  A_1 A_2 A_3;
          A2 A1  A0  A_1 A_2;
          A3 A2  A1  A0  A_1;
          A4 A3  A2  A1  A0];

B_t =    [B0 B_1 B_2 B_3 B_4;
          B1 B0  B_1 B_2 B_3;
          B2 B1  B0  B_1 B_2;
          B3 B2  B1  B0  B_1;
          B4 B3  B2  B1  B0];
   
% load('A_toep_i.mat')
% load('B_toep_i.mat')


N = (2*pi/T_x)*blkdiag(-2i*eye(4),-1i*eye(4),0i*eye(4),1i*eye(4),2i*eye(4));

% K =  [29.0282   -6.8217   -0.3859   -0.2383];
% K = [28.9452   -6.3105   -0.4698    0.0283];
% K = [71.5075  -10.8148    6.0765    2.2349];

% Optimization
% K0 = [-2.9415   16.9707    0.3163   -8.6888];
% 
% options = optimoptions(@fminunc,'MaxFunctionEvaluations',800);
% K = fminunc(@(K) objectivefcn1(K,T_angle,A_t,B_t),K0,options)
%  K = state_feedback4(T_angle,c,phase,T_x,a,i,k)

% K = [73.2581  -27.7655   35.2485   -0.6214];
K=10^3*[1.3643    0.0186    0.0106   -0.0022];
% K = 10^3*[1.9701   -0.0253    0.3274   -0.0042];
% K=10^3*[1.5    0.019    0.011   -0.003];
K=10^3*[1.5    0.019    0.017   -0.004];

% K =10^06*[2.0385   -0.0004    0.0008    0.0000];

% p=[-113.32 -1.94+1.84i -1.94-1.84i -0.25];
% p=[-113.32 -1.94+1.84i -1.94-1.84i -0.25];
% K=place(A0,B0,p)
% Q=diag([1 100 0 0]);
% [K,~,~] = lqr(A0,B0,Q,0.1,0)
% K = [15.1984   82.2251   22.3265   86.1642];
% K = [-3.5083  217.9166  -35.1181    2.9816]
% K = pso_var(T_angle,A_t,B_t);

Kh = blkdiag(K,K,K,K,K);
% 
% Acont = A_toep_i - N - B_toep_i*Kh;
% 
% eigenvalues = eig(Acont)
% eigenA0 = eig(A0-B0*K)

% F = 100;
%     load('A_t.mat')
%     load('B_t.mat')
% 
% K = [92.2150  -15.2291    6.3378    2.2090];
% K = [31.8375  -11.1479   -0.3112   -0.2723];
% K =[71.507500 -10.814800 6.076500 2.234900];

% Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;
eigenval = eig(A_t-N)

eigenvalues = eig(Acont)
figure
plot(real(eigenval),imag(eigenval),'om')
hold on
plot(real(eigenvalues),imag(eigenvalues),'*c')
hold on
plot(real(eig(A0-B0*K)),imag(eig(A0-B0*K)),'*b')
legend('açýk devre','kapalý devre','temel harmonik')
grid on

% K = [ 0 0 0 0 ];
% for j = 1:100
% %     73.2581  -27.7655   35.2485   -0.6214
%     
%     cost(j)= objectivefcn1(K,T_x,A_t,B_t);
%     K = K + [1 -0.4 0.5 -0.01];
% end
% j=1:100;
% figure
% plot(j,cost)
% grid on
% 

C=eye(20,20);
D=zeros(20,5);
sys=ss(A_t-N,B_t,C,D);
htfs=idtf(sys);

