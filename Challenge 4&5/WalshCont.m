clear;
load("test_signal1");
load("test_signal2");
step = 0.25;
k = 4
resFirst = 2^k;

U = (1./sqrt(resFirst)).*walsh(resFirst);
for i = t_start:step:t_end-step
    d = linspace(t_start,t_end,resFirst)
    c = zeros(1,resFirst);
    for s = 1:resFirst
       c(1,s)=trapz() 
        
    end
    
    
end