% variable state size PSO
function K = pso_sim(T_alpha,ca,phase,T_x,cx,i,n)
swarm_size = 32;
maxiter = 100;
fbest = 3;
fb = zeros(1,swarm_size);
c1 = 2;
c2 = 2;
n = 4;
for k = 1:n
    Kgbest(1,k) = 0;
end

for k = 1:swarm_size
    
    K(k,:) = [100*rand 100*rand-50 10*rand-5 10*rand-5];
%     K(k,:) = state_feedback(T_alpha,ca,phase,T_x,cx,i,n,K(k,:))
    v(k,:) = randn(1,n);
    fb(k) = pend_sim(K(k,:));
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
        f(k) = pend_sim(K(k,:));
        if f(k)<fb(k)
            fb(k)= f(k);
            Kbest(k,:) = K(k,:);
            if f(k)<fbest
                fbest(s) = f(k);
                Kgbest = K(k,:);
                iter(s) = k;
                s = s+1;
            end
        end
    end

end
         K = Kgbest;
         K = fminunc(@pend_sim,Kgbest);

end
