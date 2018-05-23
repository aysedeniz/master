function f = objectivefcn1(K,T_angle,A_t,B_t)

%     load('A_t.mat')
%     load('B_t.mat')

%     f1 = 0;
%     f2 = 0;
%     T_x = 1;
%     
%     N = (2*pi/T_x)*blkdiag(-2i*eye(4),-1i*eye(4),0i*eye(4),1i*eye(4),2i*eye(4));
% 
%     Kh = blkdiag(K,K,K,K,K);
%     
%     Acont = A_t - N - B_t*Kh;
%     eigen = real(eig(Acont));
%     f = max(eigen);
%     for n = 1:20
%         if (eigen(n)>-5)
%             f1 = f1 + eigen(n);
%         else 
%             f2 = f2 + eigen(n);
%         end
%     end
%     f = f1*100 + f2;

%% Inverted obj fcn
n = length(A_t)/5;
N = (2*pi/T_angle)*blkdiag(-2i*eye(n),-1i*eye(n),0i*eye(n),1i*eye(n),2i*eye(n));
Kh = blkdiag(K,K,K,K,K);

Acont = A_t - N - B_t*Kh;

eigenvalues = eig(Acont);

eigmax = max(real(eigenvalues));
eigmin = min(real(eigenvalues));
f = 4000*eigmax+0.5*eigmin+5*norm(K(3:4))+norm(K(1:2));

% if eigmax <-0.25 && K(1)>0
%     f = 0;
% else
%     f = 10+abs(eigmax+1)^2;
% end

end