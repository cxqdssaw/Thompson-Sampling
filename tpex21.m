clear
clc

N = 1;% we may choose different d for N times
T = 100;% number of episodes
H = 99;% number of timesteps
d = random('Normal',0,1,N,1);% true value of d
K = 1000;
for i = 1:N
    reward = max(d(i), 0.5);
    regret = zeros(T,K);
    regb = zeros(T,K);
    sum = 0;
    for j = 1:K
        for epi = 1:T
            rt = random('Binomial',1,0.5);
            if rt>0.5
                regret(epi,j) = reward - d(i);
                break;
            else
                regret(epi,j) = reward - 0.5;
            end
        end
        sum  = sum + epi;
    end
    fprintf('(a)learn the optimal after %f episodes\n', sum/K);
  %  plot(1:T, regret, 'b-');
    plot(1:T, mean(regret,2), 'b-');
    hold on
    sum = 0;
    for j = 1:K
       
        for epi = 1:T
            sample = random('Binomial',H,0.5);
            if sample == H
                regb(epi,j) = reward - d(i);
                break;
            else
                regb(epi,j) = reward - 0.5;
            end
        end
        sum = sum + epi;
    end
    epi = sum / K;
    if epi < T
        fprintf('(b)learn the optimal after %f episodes\n', epi);
    else
        fprintf('(b)Can NOT learn the optimal after %f episodes\n', epi);
    end
    plot(1:T, mean(regb,2), 'r-');
    %}
    
end
     
