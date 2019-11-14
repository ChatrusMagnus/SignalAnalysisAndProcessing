clear;
load("test_signal1");
load("test_signal2");
tic;
t_end = 7.5;
tic 
freq = 10000;
k = 13
resFirst = 2^k;
CoefSaving = 1;
U = (1./sqrt(resFirst)).*walsh(resFirst);
step = 0.25;
for i = t_start:step:t_end-step
    basis = linspace(i,i+step,resFirst);
    x = test_signal2(basis);
    c = U'*x;
    c(resFirst*CoefSaving:resFirst,1) = 0;
    y = U*c;
    plot(basis,c,'d');
    hold on;
    if i == t_start
        result = y;
    else
        tmp = result;
        result = cat(1,tmp,y);
    end 
end

j = linspace(t_start,t_end,length(result));
figure
plot(j,test_signal2(j)','b');
hold on;
plot(j,result,'r');
hold on;
testOri = linspace(t_start,t_end,length(result));

toc;


sound(test_signal2(testOri),length(result)./7.5);
pause(10);
sound(result,length(result)./7.5);
ene = trapz(j,result.^2);
eneOr = trapz(j,test_signal2(j)'.^2)
disp((eneOr-ene)./eneOr.*100);
