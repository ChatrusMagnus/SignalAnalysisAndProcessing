clear;
tic;
load("test_signal1");
load("test_signal2");
base = linspace(0,7.5,100000);
x = test_signal1(base);                   % Single ecg wave
y = fwht(x);                      % Fast Walsh-Hadamard transform
figure
subplot(2,1,1)
plot(base,x)
xlabel('Sample index')
ylabel('Amplitude')
title('ECG Signal')
subplot(2,1,2)
plot(abs(y))
xlabel('Sequency index')
ylabel('Magnitude')
title('WHT Coefficients')
y(30000:length(x)) = 0;            % Zeroing out the higher coefficients    
xHat = ifwht(y);                  % Signal reconstruction using inverse WHT  
figure
plot(base,x)
hold on
llll = linspace(t_start,t_end+2.344,131072);
plot(llll,xHat,'r')
xlabel('Sample index')
ylabel('ECG signal amplitude')
legend('Original Signal','Reconstructed Signal')
toc;
eneOr = integral(@(x) (test_signal1(x)').^2,t_start,7)
ene = trapz(llll,xHat.^2);
disp((eneOr-ene)./eneOr);
sound(xHat,131072./7.5)
