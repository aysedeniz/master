function [c,ceq]=nlc(K,T_angle,A_t,B_t)
n = length(A_t)/5;
N = (2*pi/T_angle)*blkdiag(-2i*eye(n),-1i*eye(n),0i*eye(n),1i*eye(n),2i*eye(n));
Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;

eigenvalues = eig(Acont);

c = max(real(eigenvalues))+0.01;
ceq =[];
