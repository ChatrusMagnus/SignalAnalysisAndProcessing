clear;
load("test_signal1");
load("test_signal2");
NOU = 10;
ene = 0;
res = 100;
a = 4.^1;
COE = zeros(1,NOU);
step = ((1./a));
lop = ((1./a)*2);
max = 7/lop;
result = zeros(1,res.*(7./(step.*2)));

tic
parfor i = 1:1:NOU
    f = @(x) legendreP(i-1,x.*a).^2;
    COE(1,i) = 1./sqrt(integral(f,-(1./a),(1./a)));
end
toc
tic
for i = 0:max
    j = i*lop;
    disp(j);
    q = linspace(j,j+(step.*2),res);
    tmpEne = zeros(NOU,res);
    for k = 1:1:NOU
        p = @(x) COE(1,k).*legendreP(k-1,(x-step-j).*a).*test_signal1(x)';
        tmp = trapz(q,p(q));
        w =@(x) tmp.*COE(1,k).*legendreP(k-1,(x-step-j).*a);  
        tmpEne(k,:) = w(q);
   end
    finalEne = ones(1,NOU)* tmpEne;
    result(1,(j.*a./2*res)+1:1:((j.*a./2 + 1)*res)) = finalEne;
    ene = ene + trapz(q,finalEne.^2);
    
    
    
end

length(result)
sl =linspace(0,7,length(result));
eneOr = integral(@(x) (test_signal1(x)').^2,t_start,7);
plot(sl,test_signal1(sl));
hold on;
plot(sl,result);
hold on;
disp((eneOr-ene).*100);
toc