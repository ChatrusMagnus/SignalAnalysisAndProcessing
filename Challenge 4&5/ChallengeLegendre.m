clear;
load("test_signal1");
load("test_signal2");
NOU = 100;
ene = 0;
res = 1000;
COE = zeros(1,NOU);
for i = 1:1:NOU
    f = @(x) legendreP(i-1,x).^2;
    COE(1,i) = 1./sqrt(integral(f,-1,1));
end
tic
for j = 0:2:4
    disp(j);
    q = linspace(j,j+2,res);
    tmpEne = zeros(NOU,res);
    for k = 1:1:NOU
        p = @(x) COE(1,k).*legendreP(k-1,x-1-j).*test_signal1(x)';  
        tmp = trapz(q,p(q));
        w =@(x) tmp.*COE(1,k).*legendreP(k-1,x-1-j);  
        tmpEne(k,:) = w(q);
   end
    finalEne = ones(1,NOU)* tmpEne ;
        plot(q,finalEne);
        grid on;
        axis([0 7.5 -1 1]);
        hold on;
    ene = ene + trapz(q,finalEne.^2);
end
plot(q,test_signal1(q));
        grid on;
        axis([0 7.5 -1 1]);
        hold on;


eneOr = integral(@(x) (test_signal1(x)').^2,t_start,6);


disp((eneOr-ene).*100);
toc