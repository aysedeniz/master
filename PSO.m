% PSO

clc
close all

swarm_size = 32;
maxiter = 100;
fbest = 0;
Kgbest = [0 0 0 0];
fb = zeros(1,swarm_size);
c1 = 2;
c2 = 2;

for k = 1:swarm_size
    
    K(k,:) = randn(1,4);
    K(k,:) = fminunc(@objectivefcn1,K(k,:));
    v(k,:) = randn(1,4);
    fb(k) = objectivefcn1(K(k,:));
    Kbest(k,:) = K(k,:);
    if fb(k)<fbest
        fbest = fb(k);
        Kgbest = K(k,:);
    end
end
s=1;
for j = 1:maxiter
    for k = 1:swarm_size
        v(k,:) = v(k,:) + c1 * rand(1,4) .* (Kbest(k,:) - K(k,:)) + c2 * rand(1,4) .* (Kgbest - K(k,:));
        K(k,:) = K(k,:) + v(k,:);
        f(k) = objectivefcn1(K(k,:));
        if f(k)<fb(k)
            fb(k)= f(k);
            Kbest(k,:) = K(k,:)
            if f(k)<fbest
                fbest(s) = f(k);
                Kgbest = K(k,:);
                iter(s) = k
                s = s+1
            end
        end
    end

end

         K = fminunc(@objectivefcn1,Kgbest)

Kgbest
plot(fbest)