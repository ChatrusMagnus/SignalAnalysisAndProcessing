 %%% This program to emulate a data transmission signal
close all;
clear all;

% number of symbols (a power of two)
Ns=2^10;

% number of time-samples (of "points") per symbol (a power of two)
Nps=2^8;

% creating the symbols
r=rand(1,Ns);
% polarity=1 --> unipolar
% polarity=2 --> bipolar
polarity=2;

if polarity==1
b=(r>0.5);
end;
if polarity==2
b=(r>0.5)-(r<0.5);
end;

% the elementary pulse used
% roll-off
alpha=0.7;
% symbol time
tau=1e-6;
% generating a single pulse
t0=tau/Nps;
tp=-tau:t0:tau;
q=root_raised_cosine(tp,alpha,tau/2,0,0);

figure;
xlabel('time','fontsize',16);
title('elementary pulse used, $q(t)$','fontsize',16,'Interpreter','LaTex');
grid on;
axis([-tau, tau, -0.1,1.1])
hold on;
plot([-tau, tau], [0 , 0],'k', 'linewidth', 1)
plot([0, 0], [0 , 1.1],'k', 'linewidth', 1)
plot(tp,q,'r', 'linewidth',2);
pause;

%% single pulse spectrum
Npf=2^8;
f0=1/tau/Npf;
f=[-5/tau:f0:5/tau];
P=tau*sinc(tau*f).*cos(pi*alpha*tau*f)./(1-(2*alpha*tau*f).^2);

figure;
xlabel('frequency','fontsize',16);
title('elementary pulse Fourier transform $Q(f)$','fontsize',16, 'Interpreter','LaTex');
grid on;
axis([-5/tau, 5/tau, -0.5*tau,1.5*tau])
hold on;
plot([-5/tau, 5/tau], [0 , 0],'k', 'linewidth', 1)
plot([0, 0], [0 , 1.5*tau],'k', 'linewidth', 1)
plot(f,P,'r','linewidth', 2);
pause;

%% generating the time-domain signal
tx_signal=zeros(1,Ns*Nps+Nps+1);
for ns=1:Ns
    nd=(ns-1)*Nps;
    tx_signal(1+nd:2*Nps+1+nd)=tx_signal(1+nd:2*Nps+1+nd)+q*b(ns);
end;

figure;
xlabel('time','fontsize',16);
title('overall transmitted signal $s(t)$','fontsize',16,'Interpreter','LaTex');
grid on;
axis([0,t0*Nps*(Ns+1), min(min(tx_signal)*1.1,-0.1),max(tx_signal)*1.1])
hold on;
plot(t0*[0:Nps*(Ns+1)],tx_signal,'r');
plot([0, t0*Nps*(Ns+1)], [0 , 0],'k', 'linewidth', 1.5)
plot([0, 0], [0 , max(tx_signal)*1.1],'k', 'linewidth', 1.5)
pause;

%% generating the Tx signal Fourier transform
tx_signal_spectrum=zeros(1,10*Npf+1);
data_spectrum=zeros(1,10*Npf+1);
for ns=1:Ns
%tx_signal_spectrum=tau*sinc(tau*f).*cos(pi*alpha*tau*f)./(1-(2*alpha*tau*f).^2).*exp(-j*2*pi*(ns*Nps)*t0*f)+tx_signal_spectrum;
tx_signal_spectrum=P.*exp(-j*2*pi*(ns*Nps)*t0*f)*b(ns)+tx_signal_spectrum;
% the raw data sequence spectrum, scaled so that it has the same height at
% the origin max(P) as the pulse transform P.
data_spectrum=max(P)*exp(-j*2*pi*(ns*Nps)*t0*f)*b(ns)+data_spectrum;

end;

figure;
xlabel('frequency','fontsize',16);
title('Abs. value of Tx signal Fourier transforms, $|D(f)|$, $|Q(f)|$ and $|S(f)|$','fontsize',16,'Interpreter','LaTex');
grid on;
hold on;
axis([-5/tau, 5/tau, -0.2*max(abs(tx_signal_spectrum)),1.1*max(abs(tx_signal_spectrum))]);
plot(f,abs(data_spectrum),'c');
plot([-5/tau, 5/tau], [0 , 0],'k', 'linewidth', 2)
plot([0, 0], [0 , max(abs(tx_signal_spectrum))],'k', 'linewidth', 2);pause;
plot(f,abs((polarity/2)*sqrt(Ns)*P),'g','Linewidth',2);pause;
plot(f,abs(tx_signal_spectrum),'r');pause;
% superimposing the pulse spectrum


figure;
xlabel('frequency','fontsize',16);
title('Tx signal Fourier transform $S(f)$','fontsize',16,'Interpreter','LaTex');
grid on;
hold on;
axis([-5/tau, 5/tau, -0.2*max(abs(tx_signal_spectrum)),1.1*max(abs(tx_signal_spectrum))]);
plot(f,abs(tx_signal_spectrum),'r')
plot([-5/tau, 5/tau], [0 , 0],'k', 'linewidth', 2)
plot([0, 0], [0 , max(abs(tx_signal_spectrum))],'k', 'linewidth', 2);pause;
% superimposing the pulse spectrum

%% generating a translated version
% wide frequency array
fw=[-50/tau:f0:50/tau];
% channel upconversion
f_up=20e6;
P_up=0.5*tau*sinc(tau*(fw-f_up)).*cos(pi*alpha*tau*(fw-f_up))./(1-(2*alpha*tau*(fw-f_up)).^2)+...
    0.5*tau*sinc(tau*(fw+f_up)).*cos(pi*alpha*tau*(fw+f_up))./(1-(2*alpha*tau*(fw+f_up)).^2);
    
tx_signal_spectrum_up=zeros(1,100*Npf+1);
for ns=1:Ns
%tx_signal_spectrum=tau*sinc(tau*f).*cos(pi*alpha*tau*f)./(1-(2*alpha*tau*f).^2).*exp(-j*2*pi*(ns*Nps)*t0*f)+tx_signal_spectrum;
tx_signal_spectrum_up=P_up.*exp(-j*2*pi*(ns*Nps)*t0*fw)*b(ns)+tx_signal_spectrum_up;

end;

h=figure;
xlabel('frequency','fontsize',16);
title('upconverted Tx signal Fourier transform','fontsize',16,'Interpreter','LaTex');
grid on;
hold on;
plot([-50/tau, 50/tau], [0 , 0],'k', 'linewidth', 2)
plot([0, 0], [0 , max(abs(tx_signal_spectrum_up))],'k', 'linewidth', 2)
plot(fw,abs(tx_signal_spectrum_up),'r');
% superimposing the signal spectrum
plot(fw,abs((polarity/2)*sqrt(Ns)*P_up),'g');
axis([-50/tau, 50/tau, -0.2*max(abs(tx_signal_spectrum_up)),1.1*max(abs(tx_signal_spectrum_up))])
pause;

%% generating a second translated version
% channel upconversion
f_up=28e6;
P_up=0.5*tau*sinc(tau*(fw-f_up)).*cos(pi*alpha*tau*(fw-f_up))./(1-(2*alpha*tau*(fw-f_up)).^2)+...
    0.5*tau*sinc(tau*(fw+f_up)).*cos(pi*alpha*tau*(fw+f_up))./(1-(2*alpha*tau*(fw+f_up)).^2);
    
tx_signal_spectrum_up=zeros(1,100*Npf+1);
for ns=1:Ns
%tx_signal_spectrum=tau*sinc(tau*f).*cos(pi*alpha*tau*f)./(1-(2*alpha*tau*f).^2).*exp(-j*2*pi*(ns*Nps)*t0*f)+tx_signal_spectrum;
tx_signal_spectrum_up=P_up.*exp(-j*2*pi*(ns*Nps)*t0*fw)*b(ns)+tx_signal_spectrum_up;

end;
plot(fw,abs(tx_signal_spectrum_up),'r');
% superimposing the signal spectrum
plot(fw,abs((polarity/2)*sqrt(Ns)*P_up),'g');
axis([-50/tau, 50/tau, -0.2*max(abs(tx_signal_spectrum_up)),1.1*max(abs(tx_signal_spectrum_up))])
figure(h);
pause;

%% generating a third translated version
% channel upconversion
f_up=40e6;
P_up=0.5*tau*sinc(tau*(fw-f_up)).*cos(pi*alpha*tau*(fw-f_up))./(1-(2*alpha*tau*(fw-f_up)).^2)+...
    0.5*tau*sinc(tau*(fw+f_up)).*cos(pi*alpha*tau*(fw+f_up))./(1-(2*alpha*tau*(fw+f_up)).^2);
    
tx_signal_spectrum_up=zeros(1,100*Npf+1);
for ns=1:Ns
%tx_signal_spectrum=tau*sinc(tau*f).*cos(pi*alpha*tau*f)./(1-(2*alpha*tau*f).^2).*exp(-j*2*pi*(ns*Nps)*t0*f)+tx_signal_spectrum;
tx_signal_spectrum_up=P_up.*exp(-j*2*pi*(ns*Nps)*t0*fw)*b(ns)+tx_signal_spectrum_up;

end;
plot(fw,abs(tx_signal_spectrum_up),'r');
% superimposing the signal spectrum
plot(fw,abs((polarity/2)*sqrt(Ns)*P_up),'g');
axis([-50/tau, 50/tau, -0.2*max(abs(tx_signal_spectrum_up)),1.1*max(abs(tx_signal_spectrum_up))])
figure(h);

