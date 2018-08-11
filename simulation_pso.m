clear
close all
clc
fl = 1; % cost function choice 1: nl sys error 2: nl-lin 3: smoothed nl-lin
swarm_size = 64;
maxiter = 100;
fb = zeros(1,swarm_size);
c1 = 2;
c2 = 2;
n = 4;
% for k = 1:n
%     Kgbest(1,k) = 0;
% end
Kgbest=[69.3726  -34.4590   20.0000    0.5004];
    if fl == 1
        [fbest,~,~] = pend_sim(Kgbest);
    elseif fl == 2
        [~,fbest,~] = pend_sim(Kgbest);
    else
        [~,~,fbest] = pend_sim(Kgbest);
    end
for k = 1:swarm_size
    
    K(k,:) = [100*rand 100*rand-50 10*rand-5 10*rand-5];
%     K(k,:) = state_feedback(T_alpha,ca,phase,T_x,cx,i,n,K(k,:))
    v(k,:) = randn(1,n);
    if fl == 1
        [fb(k),~,~] = pend_sim(K(k,:));
    elseif fl == 2
        [~,fb(k),~] = pend_sim(K(k,:));
    else
        [~,~,fb(k)] = pend_sim(K(k,:));
    end
    Kbest(k,:) = K(k,:);
    if fb(k)<fbest
        fbest = fb(k);
        Kgbest = K(k,:);
    end
end
s=1;
for j = 1:maxiter
    for k = 1:swarm_size
        v(k,:) = v(k,:) + c1 * rand(1,n) .* (Kbest(k,:) - K(k,:)) + c2 * rand(1,n) .* (Kgbest - K(k,:));
        K(k,:) = K(k,:) + v(k,:);
        if fl == 1
            [f(k),~,~] = pend_sim(K(k,:));
        elseif fl == 2
            [~,f(k),~] = pend_sim(K(k,:));
        else
            [~,~,f(k)] = pend_sim(K(k,:));
        end
        if f(k)<fb(k)
            fb(k)= f(k);
            Kbest(k,:) = K(k,:);
            if f(k)<fbest
                fbest = f(k);
                Kgbest = K(k,:);
                iter(s) = k;
                s = s+1;
            end
        end
    end
j
display('---------------------------------------------------------------------------------------------')
end
         K = Kgbest;
         display("-------------------------------------end----------------------------------------")
         K = fminunc(@pend_sim,Kgbest);
