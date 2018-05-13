%% fourier decomposition
function K = state_feedback4(T_alpha,ca,phase,T_x,cx,i,n)

T = T_x;
Fs = 400;
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
% c=(0:(Ls-1))/Fs;
% figure
% plot(c,alpha)
% figure
% plot(c,alpha_d)
% figure
% plot(c,x_d)
% figure
% plot(c,V)
% figure(90)
% plot((0:(Ls-1))/Fs,A_n(1,:))
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

K = [73.2581  -27.7655   35.2485   -0.6214];

Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;
% eigenval = eig(A_t)
% Optimization

% K = pso_var(T,A_t,B_t);
% K = fminunc(@(K) objectivefcn1(K,T,A_t,B_t),[73.2581  -27.7655   35.2485   -0.6214]);
Q=diag([1 100 0 0]);
[K,~,~] = lqr(A0,B0,Q,0.001,0);
Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;
eigenval = eig(A_t-N)

eigenvalues = eig(Acont)
% figure
% plot(real(eigenval),imag(eigenval),'om','MarkerSize',10)
% hold on
% plot(real(eigenvalues),imag(eigenvalues),'*c','MarkerSize',10)
% % hold on
% % plot(real(eig(A0)),imag(eig(A0)),'*b')
% legend('açýk devre','kapalý devre','temel harmonik')
% grid on

% grid on
end
