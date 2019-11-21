function Pf=root_raised_cosine(f,beta,Rs,center_freq,root)

% f is the array of frequencies at which the function is evaluated
% beta is the roll-off factor
% Rs is the Baud rate
% center_freq is the channel center frequency
% root at 1 indicates that a root-raised-cosine spectrum must be returned
% root at 0 indicates that a rasied-cosine spectrum must be returned

Ts=Rs^-1;
passband=(1-beta)/2/Ts;
stopband=(1+beta)/2/Ts;

% initializing
Pf=zeros(1,length(f));

for nch=1:length(center_freq)
    
    ff=abs(f-center_freq(nch));
    tf=abs(f-center_freq(nch))-passband;

if beta==0
    Pf=(tf<=0)+Pf;
else
    if root==1
    Pf=...
    (tf<=0)+...
    sqrt(1/2*(1+cos(pi*Ts/beta*(tf))).*(tf>0).*(abs(ff)<=stopband))+Pf;
    else
    Pf=...
    (tf<=0)+...
    (1/2*(1+cos(pi*Ts/beta*(tf))).*(tf>0).*(abs(ff)<=stopband))+Pf;
    end;
end;

end;

return





return
t=[-3*pi:0.05:3*pi]
alpha=0.
alpha=0:0.2:1
for i=1:6
T=1;
rcos(i,:)=sin(pi*t/T)./(pi*t/T).*cos(pi*t/T*alpha(i))./(1-4*alpha(i)^2/T^2.*t.^2);
end;
plot(t,rcos,'linewidth',2.5);
grid on;


