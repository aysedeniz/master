% Frequency search
close all
clc
clear

for j = 1:20
    
    K = [69.3726  -34.4590   20.0000    0.5004];
    T = j*0.02;
    [err(j), conv_err(j)] = pend_sim(K,T);

k = 1;
n = 1;
T_alpha = T;
ca = 0.1;
[amp, T_x, phase, ofs] = find_position(ca, T_alpha, k);
cx = amp;
phase_ofs = 0;

Fs = 72/T;
Ls = Fs*T*100;

dx = 2*pi/T_x;
dalpha = 2*pi/T_alpha;

c = 0;
for t = (0:(Ls-1))/Fs  
   c = c+1;
   alpha(c) = pi*n+ca*sin(dalpha*t+phase*abs(i-1));
   alpha_d(c) = ca*dalpha*cos(dalpha*t+phase*abs(i-1));
   x_d(c) = cx*dx*cos(dx*t+phase*i);
   V(c) = openloop_input(cx, T_x, ca, T_alpha,i,n, phase,0, t);
   [A,B]=linearized_model(alpha(c),x_d(c),alpha_d(c),V(c));
    for k=0:1
        for j=1:3
            A_n(j+(k*3),c) = A(3+k,j+1);
            if j==1
                B_n(k+1,c) = B(k+3);
            end
        end
    end
end

for k=1:6
    f(k,:) = fft(A_n(k,:))/Ls;
    f1(k,:) = 2*abs(f(k,:));
    f1(k,1) = f1(k,1)/2;
    f2(k,:) = f1(k,1:Ls/2+1);
end


for k=1:2
    fb(k,:) = fft(B_n(k,:))/Ls;
    fb1(k,:) = 2*abs(fb(k,:));
    fb1(k,1) = fb1(k,1)/2;
    fb2(k,:) = fb1(k,1:Ls/2+1);
end

loc(1) = (Ls/(T*Fs))+1;
loc(2) = 2*(Ls/(T*Fs))+1;
loc(3) = 3*(Ls/(T*Fs))+1;
loc(4) = 4*(Ls/(T*Fs))+1;
loc(5) = 5*(Ls/(T*Fs))+1;

for k=1:6
    a1(k) = f2(k,loc(1))*exp(angle(f(k,loc(1)))*1i)/2;
    a_1(k) = f2(k,loc(1))*exp(-angle(f(k,loc(1)))*1i)/2;
    a2(k) = f2(k,loc(2))*exp(angle(f(k,loc(2)))*1i)/2;
    a_2(k) = f2(k,loc(2))*exp(-angle(f(k,loc(2)))*1i)/2;
    a3(k) = f2(k,loc(3))*exp(angle(f(k,loc(3)))*1i)/2;
    a_3(k) = f2(k,loc(3))*exp(-angle(f(k,loc(3)))*1i)/2;
    a4(k) = f2(k,loc(4))*exp(angle(f(k,loc(4)))*1i)/2;
    a_4(k) = f2(k,loc(4))*exp(-angle(f(k,loc(4)))*1i)/2;
end


A0 = [0      0      1      0
      0      0      0      1;
      0      f(1,1) f(2,1) f(3,1);
      0      f(4,1) f(5,1) f(6,1)];
  
B0 = [0; 0; fb(1,1); fb(2,1)];

A1 = [0     0       0      0
      0     0       0      0;
      0     a1(1)   a1(2)  a1(3);
      0     a1(4)   a1(5)  a1(6)];
  
A_1 = [0     0       0      0
       0     0       0      0;
       0     a_1(1)  a_1(2) a_1(3);
       0     a_1(4)  a_1(5) a_1(6)];
  
B1 = [0; 0; fb2(1,loc(1))*exp(angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(angle(fb(2,loc(1)))*1i)/2];
B_1 = [0; 0; fb2(1,loc(1))*exp(-angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(-angle(fb(2,loc(1)))*1i)/2];


A2 = [0     0       0      0
      0     0       0      0;
      0     a2(1)   a2(2)  a2(3);
      0     a2(4)   a2(5)  a2(6)];
  
A_2 = [0     0       0      0
       0     0       0      0;
       0     a_2(1)  a_2(2) a_2(3);
       0     a_2(4)  a_2(5) a_2(6)];

B2 = [0; 0; fb2(1,loc(2))*exp(angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(angle(fb(2,loc(2)))*1i)/2];
B_2 = [0; 0; fb2(1,loc(2))*exp(-angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(-angle(fb(2,loc(2)))*1i)/2];

A3 = [0     0       0      0
      0     0       0      0;
      0     a3(1)   a3(2)  a3(3);
      0     a3(4)   a3(5)  a3(6)];
  
A_3 = [0     0       0      0
       0     0       0      0;
       0     a_3(1)  a_3(2) a_3(3);
       0     a_3(4)  a_3(5) a_3(6)];

B3 = [0; 0; fb2(1,loc(3))*exp(angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(angle(fb(2,loc(3)))*1i)/2];
B_3 = [0; 0; fb2(1,loc(3))*exp(-angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(-angle(fb(2,loc(3)))*1i)/2];

A4 = [0     0       0      0
      0     0       0      0;
      0     a4(1)   a4(2)  a4(3);
      0     a4(4)   a4(5)  a4(6)];
  
A_4 = [0     0       0      0
       0     0       0      0;
       0     a_4(1)  a_4(2) a_4(3);
       0     a_4(4)  a_4(5) a_4(6)];

B4 = [0; 0; fb2(1,loc(4))*exp(angle(fb(1,loc(4)))*1i)/2;...
    fb2(2,loc(4))*exp(angle(fb(2,loc(4)))*1i)/2];
B_4 = [0; 0; fb2(1,loc(4))*exp(-angle(fb(1,loc(4)))*1i)/2;...
    fb2(2,loc(4))*exp(-angle(fb(2,loc(4)))*1i)/2];

    
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
    
N = (2*pi/T)*blkdiag(-2i*eye(4),-1i*eye(4),0i*eye(4),1i*eye(4),2i*eye(4));
Kh = blkdiag(K,K,K,K,K);
Acont = A_t - N - B_t*Kh;
eigenval = eig(A_t-N);

eigenvalues = eig(Acont);

figure
plot(real(eigenval),imag(eigenval),'om','MarkerSize',10)
hold on
plot(real(eigenvalues),imag(eigenvalues),'*c','MarkerSize',10)
hold on
plot(real(eig(A0)),imag(eig(A0)),'*b')
legend('a�?k devre','kapal? devre','temel harmonik')
grid on
end

figure
subplot(2,1,1)
plot(err)
grid on
subplot(2,1,2)
plot(conv_err)
grid on