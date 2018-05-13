%% fourier decomposition
function K = state_feedback(T_alpha,calpha,phase_alpha,T_xdot,cxdot,phase_xdot,xdot_ofs...
    ,T_alphadot,calphadot,phase_alphadot,T_v,c_v,phase_v,V_ofs)

T = round(10*max([T_alpha T_xdot T_alphadot]))/10
Fs = 50;
Ls = Fs*5;

    dxdot = 2*pi/T_xdot;
    dalpha = 2*pi/T_alpha;
    dalphadot = 2*pi/T_alphadot;
    d_v = 2*pi/T_v;
    

for k=0:1
    for j=1:3
        c = 0;
        for t = 0:1/Fs:(Ls-1)/Fs
            c = c+1;
            alpha = pi+calpha*sin(dalpha*t+phase_alpha);
            alpha_d = calphadot*sin(dalphadot*t+phase_alphadot);
            x_d = xdot_ofs+cxdot*sin(dxdot*t+phase_xdot);
            V = V_ofs+c_v*sin(d_v*t+phase_v);

            [A,B]=linearized_model_3state(alpha,x_d,alpha_d,V);  
            A_n(j+(k*3),c) = A(2+k,j)
            if j==1
                B_n(k+1,c) = B(k+2);
            end
        end
        j
    end
    k
end


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
    f2(k,:) = f1(k,1:Ls/2+1)
    
%     subplot(2,3,k)
%     plot(Fs*(0:(Ls/2))/Ls,f2(k,:))
end

%     figure(4)
for k=1:2
    fb(k,:) = fft(B_n(k,:))/Ls;
    fb1(k,:) = 2*abs(f(k,:));
    fb1(k,1) = fb1(k,1)/2;
    fb2(k,:) = fb1(k,1:Ls/2+1);
    
%     subplot(2,1,k)
%     plot(Fs*(0:(Ls/2))/Ls,fb2(k,:))
end

loc(1) = (Ls/(T*Fs))+1
loc(2) = 2*(Ls/(T*Fs))+1
loc(3) = 3*(Ls/(T*Fs))+1
loc(4) = 4*(Ls/(T*Fs))+1
loc(5) = 5*(Ls/(T*Fs))+1

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


A0 = [0       0      1;
      f(1,1) f(2,1) f(3,1);
      f(4,1) f(5,1) f(6,1)];
  
B0 = [0; fb(1,1); fb(2,1)];

A1 = [0       0      0;
      a1(1)   a1(2)  a1(3);
      a1(4)   a1(5)  a1(6)];
  
A_1 = [0       0      0;
       a_1(1)   a_1(2)  a_1(3);
       a_1(4)   a_1(5)  a_1(6)];
  
B1 = [0; fb2(1,loc(1))*exp(angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(angle(fb(2,loc(1)))*1i)/2];
B_1 = [0; fb2(1,loc(1))*exp(-angle(fb(1,loc(1)))*1i)/2;...
    fb2(2,loc(1))*exp(-angle(fb(2,loc(1)))*1i)/2];


A2 = [0       0      0;
      a2(1)   a2(2)  a2(3);
      a2(4)   a2(5)  a2(6)];
  
A_2 = [0       0       0;
      a_2(1)   a_2(2)  a_2(3);
      a_2(4)   a_2(5)  a_2(6)];

B2 = [0; fb2(1,loc(2))*exp(angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(angle(fb(2,loc(2)))*1i)/2];
B_2 = [0; fb2(1,loc(2))*exp(-angle(fb(1,loc(2)))*1i)/2;...
    fb2(2,loc(2))*exp(-angle(fb(2,loc(2)))*1i)/2];

A3 = [0       0      0;
      a3(1)   a3(2)  a3(3);
      a3(4)   a3(5)  a3(6)];
  
A_3 = [0       0       0;
       a_3(1)  a_3(2)  a_3(3);
       a_3(4)  a_3(5)  a_3(6)];

B3 = [0; fb2(1,loc(3))*exp(angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(angle(fb(2,loc(3)))*1i)/2];
B_3 = [0; fb2(1,loc(3))*exp(-angle(fb(1,loc(3)))*1i)/2;...
    fb2(2,loc(3))*exp(-angle(fb(2,loc(3)))*1i)/2];

A4 = [0       0      0;
      a4(1)   a4(2)  a4(3);
      a4(4)   a4(5)  a4(6)];
  
A_4 = [0       0       0;
       a_4(1)  a_4(2)  a_4(3);
       a_4(4)  a_4(5)  a_4(6)];

B4 = [0; fb2(1,loc(4))*exp(angle(fb(1,loc(4)))*1i)/2;...
    fb2(2,loc(4))*exp(angle(fb(2,loc(4)))*1i)/2];
B_4 = [0; fb2(1,loc(4))*exp(-angle(fb(1,loc(4)))*1i)/2;...
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
    
N = (2*pi/T)*blkdiag(-2i*eye(3),-1i*eye(3),0i*eye(3),1i*eye(3),2i*eye(3));

% Optimization

K = pso_var(T,A_t,B_t);

Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;
eigenval = eig(A_t);

eigenvalues = eig(Acont)
figure
plot(real(eigenval),imag(eigenval),'om')
hold on
plot(real(eigenvalues),imag(eigenvalues),'*c')
hold on
plot(real(eig(A0)),imag(eig(A0)),'*b')
legend('açýk devre','kapalý devre','temel harmonik')
grid on

% grid on
end
