function [c,ceq]=const(K)
load('A_toep.mat')
load('B_toep.mat')
T_x = 0.2;
N = (2*pi/T_x)*blkdiag(-2i*eye(4),-1i*eye(4),0i*eye(4),1i*eye(4),2i*eye(4));

Kh = blkdiag(K,K,K,K,K);
Acont = A_toep - N - B_toep*Kh;

eigen = real(eig(Acont));
c = min(eigen)-10*max(eigen);
ceq = [];


end
