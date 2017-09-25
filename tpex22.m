clear
clc

N = 10;
T = 200;
K = 100;
theta = random('Normal',0,1,N,N);

Mean = zeros(N);
variance = ones(N);

theta_est = random('Normal',Mean,sqrt(variance));
x = zeros(N);

epsilon = [0.1 0.2 0.3];

v = ones(N,1);
A1 = repmat(diag(v),1,N);
A2 = zeros(size(A1));
for i = 1:N
    A2(i,(i-1)*N+1:i*N) = v;
end
Aeq = [A1; A2];
beq = ones(2*N,1);
lb = zeros(N*N,1);
ub = ones(N*N,1);
intcon = 1:N*N;
options = optimoptions('linprog','Algorithm','interior-point-legacy');

x_star = linprog(-theta(:),[],[],Aeq,beq,lb,[],options);
x_star = round(x_star);
reward = sum(x_star .* theta(:));
regret = ones(T,K);

mean_greedy = Mean;
variance_greedy = variance;
theta_greedy = mean_greedy;
regret_greedy = ones(T,K);

mean_ucb = Mean;
variance_ucb = variance;
theta_ucb = mean_ucb + 2.* sqrt(variance_ucb);
regret_ucb = ones(T,K);

mean_ep = zeros(N*N,3);
variance_ep = ones(N*N,3);
theta_ep = mean_ep ;
regret_ep = ones(T,3,K);
for j = 1:K
    for i = 1:T
        
        x = linprog(-theta_est(:),[],[],Aeq,beq,lb,[],options);
        x = round(x);
        b = find(x==1);
        y = random('Normal',0,1,N,1) + theta(b);
        regret(i,j) = reward - sum(theta(b));
        Mean(b) = (Mean(b) + variance(b).* y) ./ (variance(b)+ones(N,1));%question
        variance(b) = variance(b) ./ (variance(b)+ones(N,1));
        theta_est = random('Normal',Mean,sqrt(variance));
        %}
        
        for e = 1: 3
            x_ep = intlinprog(-theta_ep(:,e),intcon,[],[],Aeq,beq,lb,[]);
            x_ep = round(x_ep);
            b_ep = find(x_ep ==1);
            y_ep = random('Normal',0,1,N,1) + theta(b_ep);
            regret_ep(i,e,j) = reward - sum(theta(b_ep));
            mean_ep(b_ep,e) = (mean_ep(b_ep,e) +...
                variance_ep(b_ep,e).* y_ep) ./ (variance_ep(b_ep,e)+ones(N,1));%question
            variance_ep(b_ep,e) = variance_ep(b_ep,e)...
                ./ (variance_ep(b_ep,e)+ones(N,1));
            flag = random('Unif',0,1);
            if flag > epsilon(e)
                theta_ep(:,e) = mean_ep(:,e);
            else
                theta_ep(:,e) = random('Normal',mean_ep(:,e),variance_ep(:,e));%sample randomly
            end
            
        end
        %}
        x_greedy = intlinprog(-theta_greedy(:),intcon,[],[],Aeq,beq,lb,[]);
        x_greedy = round(x_greedy);
        b_greedy = find(x_greedy ==1);
        y_greedy = random('Normal',0,1,N,1) + theta(b_greedy);
        regret_greedy(i,j) = reward - sum(theta(b_greedy));
        mean_greedy(b_greedy) = (mean_greedy(b_greedy) +...
            variance_greedy(b_greedy).* y_greedy) ./ (variance_greedy(b_greedy)+ones(N,1));%question
        variance_greedy(b_greedy) = variance_greedy(b_greedy)...
            ./ (variance_greedy(b_greedy)+ones(N,1));
        theta_greedy = mean_greedy;
        
        x_ucb = intlinprog(-theta_ucb(:),intcon,[],[],Aeq,beq,lb,[]);
        x_ucb = round(x_ucb);
        b_ucb = find(x_ucb ==1);
        y_ucb = random('Normal',0,1,N,1) + theta(b_ucb);
        regret_ucb(i,j) = reward - sum(theta(b_ucb));
        mean_ucb(b_ucb) = (mean_ucb(b_ucb) +...
            variance_ucb(b_ucb).* y_ucb) ./ (variance_ucb(b_ucb)+ones(N,1));%question
        variance_ucb(b_ucb) = variance_ucb(b_ucb)...
            ./ (variance_ucb(b_ucb)+ones(N,1));
        theta_ucb = mean_ucb + 2.*sqrt(variance_ucb);
  
    end
end

plot(1:T, mean(regret,2), 'b-');
hold on
plot(1:T, mean(regret_greedy, 2), 'r-');
plot(1:T, mean(regret_ucb, 2), 'k-');

plot(1:T, mean(regret_ep(:,1,:),3), 'g-',1:T, mean(regret_ep(:,2,:),3), 'm-',...
    1:T, mean(regret_ep(:,3,:),3), 'c-');
 xlabel('time');
 ylabel('E[regret]');
 legend('TS','greedy','ucb','0.1-greedy','0.2-greedy','0.3-greedy');
%}
