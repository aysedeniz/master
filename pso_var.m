% variable state size PSO
function Kfin = pso_var(T_angle,A_t,B_t,i)
swarm_size =32;
maxiter = 20;
fbest = 5;
fb = zeros(1,swarm_size);
c1 = 2;
c2 = 5;
n = length(A_t)/5;
lb = [0 -100 -20 -10];
ub = [2000 100 20 10];
for k = 1:n
    Kgbest(1,k) = 0;
end

for k = 1:swarm_size
    
    K(k,:) = [10000*randn 100*randn 100*randn 10*randn];
    K(k,:) = fmincon(@(K) objectivefcn1(K,T_angle,A_t,B_t),K(k,:),[],[],[],[],lb,ub,@(K) nlc(K,A_t,B_t,T_angle,i));
    v(k,:) = randn(1,n);
    fb(k) = objectivefcn1(K(k,:),T_angle,A_t,B_t);
    Kbest(k,:) = K(k,:);
    [c,ceq] = nlc(K(k,:),A_t,B_t,T_angle,i);
    c
    const = any(c>0)
    if fb(k)<fbest && const == 0
        fbest = fb(k);
        Kgbest = K(k,:);
    end
end
s=1;
for j = 1:maxiter
    for k = 1:swarm_size
        v(k,:) = v(k,:) + c1 * randn(1,n) .* (Kbest(k,:) - K(k,:)) + c2 * randn(1,n) .* (Kgbest - K(k,:));
        K(k,:) = K(k,:) + v(k,:);
        K(k,:) = fmincon(@(K) objectivefcn1(K,T_angle,A_t,B_t),K(k,:),[],[],[],[],lb,ub,@(K) nlc(K,A_t,B_t,T_angle,i));
        f(k) = objectivefcn1(K(k,:),T_angle,A_t,B_t);
        [c,ceq] = nlc(K(k,:),A_t,B_t,T_angle,i);
        c
        const = any(c>0)
        if f(k)<fb(k) && const == 0
            fb(k)= f(k);
            Kbest(k,:) = K(k,:);
            if f(k)<fbest
                fbest(s) = f(k);
                Kgbest = K(k,:)
                iter(s) = k;
                s = s+1;
            end
        end
    end
j
end
         Kfin = Kgbest
%          K = fminunc(@(K) objectivefcn1(K,T_angle,A_t,B_t),Kgbest);

% K = fmincon(@(K) objectivefcn1(K,T_angle,A_t,B_t),K0,[],[],[],[],[],[],@(K) nlc(K,T_angle,A_t,B_t))
end

