clear all; clc; close all

mic=2:6;
% Load squared pressures measurements
P_sans = load('P_sans1.txt'); 
freq = P_sans(:,1);
P_sans = sum(P_sans(:,mic),2);

P_sans2 = load('P_sans2.txt');
P_sans2 = sum(P_sans2(:,mic),2);


meas_D = load('P_equipe3.txt');
meas_D = sum(meas_D(:,mic),2);




%% Trapz method to calculate global IL
M=15;%box mass

IL_D=trapz(freq(7:end),P_sans2(7:end))./trapz(freq(7:end),meas_D(7:end));% Il par frequence double paroi
IL_D_dB=10*log10(IL_D);% IL par frequence dB
IL_D_dB_M=10*log10(IL_D)-20*log10(M(1));% IL double pondéré par la masse


%% Control plot niveau sonore
figure(1)
semilogx(freq,10*log10(P_sans/2e-5),'+-')
hold on
semilogx(freq,10*log10(P_sans2/2e-5),'+-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB_{spl}]')
xlim([25 18000])
legend('meas1','meas2')
title('Pression sans boite')

figure(2)
semilogx(freq,10*log10(meas_D/2e-5),'+-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB_{spl}]')
xlim([25 18000])
legend('Double','Triple')
title('Pression avec boite')
%% Plot IL omega
figure(3)
semilogx(freq,10*log10(P_sans2./meas_D),'+-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
xlim([25 18000])
legend(sprintf('Double, IL = %.1f dB, IL_M = %.1f dB',IL_D_dB,IL_D_dB_M))
title('IL(\omega)')
xlim([25 6000])

 FigurePlacecement()
