clear;
load("test_signal1");
load("test_signal2");
Noo = 10;
step = (t_end-t_start)./Noo;
l = step./2;
syms h
%normalization equation
eq = -1+(2/5)*(l)^5 + h*(-(4/3)*l^3) + (h)^2*(2*l) == 0;
tmp = double(solve(eq,h));
b = tmp(2,1);
MyPol = @(x,a) -((x-a-l).^2) + b;
ene = 0;
tic
for k =0:Noo-1
 
  left = t_start+(k.*step);
  right = t_start+((k+1).*step);
  y = linspace(left,right,100);
  f =@(x) MyPol(x,k.*step).*test_signal1(x)';
  plot(y,MyPol(y,k.*step));
  hold on;
  r =integral(f,left,right,'RelTol', 1e-3);
    if k == 0
        
        result = (MyPol(y,k.*step)).*r;
    else
        tmp = result; 
        
        result = cat(2,tmp,(MyPol(y,k.*step)).*r);
    end 
  w =@(x) ((MyPol(x,k.*step)).*r).^2;
  ene = ene + integral(w,left,right,'RelTol', 1e-3); 
end

j = linspace(t_start,t_end,length(result));
figure
plot(j,test_signal1(j)','b');
hold on;
plot(j,result,'r');
hold on;
%sound(result,length(result)./7.5);
eneOr = integral(@(x) (test_signal1(x)').^2,t_start,t_end);
disp((eneOr-ene).*100);
toc