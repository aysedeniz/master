% variable state size PSO
function K = pso_var(T_angle,A_t,B_t)
swarm_size = 64;
maxiter = 100;
fbest = 3;
fb = zeros(1,swarm_size);
c1 = 2;
c2 = 2;
n = length(A_t)/5;
for k = 1:n
    Kgbest(1,k) = 0;
end

for k = 1:swarm_size
    
    K(k,:) = 100*rand(1,n)-50;
    K(k,:) = fminunc(@(K) objectivefcn1(K,T_angle,A_t,B_t),K(k,:));
    v(k,:) = randn(1,n);
    fb(k) = objectivefcn1(K(k,:),T_angle,A_t,B_t);
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
        f(k) = objectivefcn1(K(k,:),T_angle,A_t,B_t);
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
%          K = fminunc(@(K) objectivefcn1(K,T_angle,A_t,B_t),Kgbest);

K = fmincon(@(K) objectivefcn1(K,T_angle,A_t,B_t),K,[],[],[],[],[],[],@(K) nlc(K,T_angle,A_t,B_t))
end

