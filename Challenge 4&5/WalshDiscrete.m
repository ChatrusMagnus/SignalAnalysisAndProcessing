clear;
load("test_signal1");
load("test_signal2");
tic;
t_end = 7.5;
tic 
freq = 1000;
resFirst = 2^13;
CoefSaving = 1;
basis = linspace(t_start,t_end,resFirst);
x = test_signal1(basis);
U = (1./sqrt(resFirst)).*hadamard(resFirst);
c = U'*x;
c(1,freq*t_end : 2^10) = 0;
y = U*c;

step = 0.5;
for i = t_start:step:t_end-0.5
    basis = linspace(i,i+step,resFirst);
    x = test_signal1(basis);
    U = (1./sqrt(resFirst)).*hadamard(resFirst);
    c = U'*x;
    c(resFirst*CoefSaving:resFirst,1) = 0;
    y = U*c;
    if i == t_start
        result = y;
    else
        tmp = result;
        result = cat(1,tmp,y);
    end 
end

j = linspace(t_start,t_end,length(result));
length(y);
plot(j,test_signal1(j)','b');
hold on;
plot(j,result,'r');
hold on;

toc;
sound(result,length(result)./7.5);
ene = trapz(result.^2)
eneOr = trapz(test_signal1(j)'.^2)
disp((eneOr-ene)./eneOr.*100);