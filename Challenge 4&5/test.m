NOU = 5;
COE = zeros(1,NOU);
q = linspace(-0.25,0.25,1000);
for i = 1:1:NOU
    f = @(x) legendreP(i-1,x.*4);
    plot(q,f(q));
    hold on;
end
