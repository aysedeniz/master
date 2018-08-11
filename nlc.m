function [c,ceq] = nlc(K,A_t,B_t,T_angle,i)
i
    n = length(A_t)/5;
    N = (2*pi/T_angle)*blkdiag(-2i*eye(n),-1i*eye(n),0i*eye(n),1i*eye(n),2i*eye(n));
    Kh = blkdiag(K,K,K,K,K);

    Acont = A_t - N - B_t*Kh;

    eigenvalues = eig(Acont);

    eigmax = max(real(eigenvalues));
    eigmin = min(real(eigenvalues));
    Q = [0 0 0 0; 0 0 0 0; 0 0 1 0; 0 0 0 10];
    R = [10 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0];

    ceq = [];
%     c(1) = -(eigmin+50);
if i==1
    c(1) = (K*Q*K'-700);
    c(2) = (-K*R*K'+10000000);
    c(3) = (eigmin+15);
%     c(2) = eigmax+0.2;
else
    c(1) = -(eigmin+50);
end
end