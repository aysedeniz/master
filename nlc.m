function [c,ceq] = nlc(K,A_t,B_t,T_angle)
    n = length(A_t)/5;
    N = (2*pi/T_angle)*blkdiag(-2i*eye(n),-1i*eye(n),0i*eye(n),1i*eye(n),2i*eye(n));
    Kh = blkdiag(K,K,K,K,K);

    Acont = A_t - N - B_t*Kh;

    eigenvalues = eig(Acont);

    eigmax = max(real(eigenvalues));
    eigmin = min(real(eigenvalues));
    ceq = [];
    c(1) = eigmin+80;
    c(2) = eigmax+0.4;
    c(3) = norm(K(1:3))-80;
end