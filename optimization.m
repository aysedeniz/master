close all
clear
clc

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];

%% This part is to draw surfaces of cost function
    load('A_toep.mat')
    load('B_toep.mat')
    
    K0 = [0 0 0 0];

    T_x = 0.2;
    
    N = (2*pi/T_x)*blkdiag(-2i*eye(4),-1i*eye(4),0i*eye(4),1i*eye(4),2i*eye(4));
    
    for j = -100:100
        for k = 0:100
        K = K0 + [k/10 j/10 0 0];
        Kh = blkdiag(K,K,K,K,K);
        Acont = A_toep - N - B_toep*Kh;
        eigen = real(eig(Acont));
       
        f(k+1,j+101) = max(eigen);
  
        end
    end
    
    K0 = [0 0 10 0];
    for j = -100:100
        for k = 0:100
        K = K0 + [k/10 j/10 0 0];
        Kh = blkdiag(K,K,K,K,K);
        Acont = A_toep - N - B_toep*Kh;
        eigen = real(eig(Acont));
       
        f1(k+1,j+101) = max(eigen);
  
        end
    end
    fmin = min(min(f));
    fmax = max(max(f));
    f1min = min(min(f1));
    f1max = max(max(f1));
    
    C(:,:,1) = zeros(size(f)); % red
    C(:,:,2) = (f-fmin)/(fmax-fmin); % green
    C(:,:,3) = zeros(size(f)); % blue
    C1(:,:,1) = (f1-f1min)/(f1max-f1min); % red
    C1(:,:,2) = zeros(size(f1)); % green
    C1(:,:,3) = zeros(size(f1)); % blue
    [X,Y] = meshgrid(-10:0.1:10,0:0.1:10);
    surf(X,Y,f,C,'EdgeColor','none','FaceColor','interp')
    hold on
    surf(X,Y,f1,C1,'EdgeColor','none','FaceColor','interp')
%%
    
% K0 = [30 -5 30 -20];
% options = optimoptions('fmincon','MaxIter',40,'MaxFunEvals',400,'TolFun',0.001);
% tic
% [K,fval,exitflag,output] = fmincon(@objectivefcn,K0,A,b,Aeq,beq,lb,ub,@eigenval,options);
% toc
% K = fminunc(@objectivefcn1,K0);
% K


% ans = struct with fields:
%         RiseTime: 0.4000
%     SettlingTime: 2.8000
%      SettlingMin: -0.6724
%      SettlingMax: -0.5188
%        Overshoot: 24.6476
%       Undershoot: 11.1224
%             Peak: 0.6724
%         PeakTime: 1

