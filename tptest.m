
K = 3;
N = 1000;
T = 1000;
count = zeros(T,3);
for k = 1:N
    p = random('Unif',0,1,K,1);
    p = sort(p);
    par = ones(K,2);
    %time = zeros(T,1);
    for i = 1: T
        r = random('Beta', par(:,1), par(:,2));
        [pp, idx] = max(r);
       % time(i) = idx;
        count(i,idx) = count(i,idx) + 1;
        rt = random('Binomial',1,p(idx));
        par(idx,:) = par(idx,:) + [rt, 1-rt];
    end
end
count = count ./ [sum(count(1,:)) sum(count(2,:)) sum(count(3,:))];
plot(1:T, count(:,1),'b-', 1:T, count(:,2),'g-', 1:T, count(:,3),'r-');
%{
x = 0:.05:1;
figure
hold on
for j = 1: K
    pd = makedist('Beta','a',par(j,1),'b',par(j,2));
    pdfx = pdf(pd,x);
    plot(x,pdfx);
end
hold off
figure(2)
plot(1:T, time, 'g*');
%}



%{
count = zeros(K,1);
gamma = 0.95;
flag = 1;
time = 1;
obj = 0;

while flag == 1
    [p, idx] = max(par(:,1) ./(par(:,1) + par(:,2)) );
    if p >= 0.5 
        obj = obj + gamma^time * (2*p-1);
        count(idx) = count(idx) + 1;
    end
    time = time + 1;
end

%}

