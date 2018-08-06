function [c,ceq] = nlc(K,A_t,B_t,T_angle)
    n = length(A_t)/5;
    N = (2*pi/T_angle)*blkdiag(-2i*eye(n),-1i*eye(n),0i*eye(n),1i*eye(n),2i*eye(n));
    Kh = blkdiag(K,K,K,K,K);

    Acont = A_t - N - B_t*Kh;

    eigenvalues = eig(Acont);

    eigmax = max(real(eigenvalues));
    eigmin = min(real(eigenvalues));
    Q = [0 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 10];
    ceq = [];
%     c(1) = -(eigmin+50);
    c(1) = K*Q*K'-100;
%     c(2) = eigmax+0.9;
%     c(3) = norm(K(2:3))-20;
%     c(1) = norm(K(4))-3;
%     c(3) = norm(K(4))-5;
%     c(4) = norm(K(3))-30;
%     c(5) = norm(K(2))-30;
end