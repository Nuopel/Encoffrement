clear all; clc; close all

mic=2:6;
% Load squared pressures measurements
P_sans = load('P_sans1.txt'); 
freq = P_sans(:,1);
P_sans = sum(P_sans(:,mic),2);

P_sans2 = load('P_sans2.txt');
P_sans2 = sum(P_sans2(:,mic),2);


meas_D = load('P_avec_equipe2_double.txt');
% meas_D = load('P_equipe1.txt');
meas_D = sum(meas_D(:,mic),2);


meas_T = load('P_avec_equipe2_triple.txt');
meas_T = sum(meas_T(:,mic),2);

load('TC_double.mat')
load('TC_triple.mat')

%% Trapz method to calculate global IL
M=[16.86 32.54];%box mass

IL_D=trapz(freq(8:end),P_sans2(8:end))./trapz(freq(8:end),meas_D(8:end));% Il par frequence double paroi
IL_D_dB=10*log10(IL_D);% IL par frequence dB
IL_D_dB_M=10*log10(IL_D)-20*log10(M(1));% IL double pondéré par la masse

IL_T=trapz(freq(8:end),P_sans2(8:end))/trapz(freq(8:end),meas_T(8:end));
IL_T_dB=10*log10(IL_T);% IL par frequence dB triple paroi dB
IL_T_dB_M=10*log10(IL_T)-20*log10(M(2)); % IL triple pondéré par la masse

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
hold on
semilogx(freq,10*log10(meas_T/2e-5),'+-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB_{spl}]')
xlim([25 18000])
legend('Double','Triple')
title('Pression avec boite')
%% Plot IL omega
figure(3)
semilogx(freq,10*log10(P_sans2./meas_D),'+-')
hold on
semilogx(freq,10*log10(P_sans2./meas_T),'+-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
xlim([25 18000])
legend(sprintf('Double, IL = %.1f dB, IL_M = %.1f dB',IL_D_dB,IL_D_dB_M),sprintf('Triple, IL = %.1f dB, IL_M = %.1f dB',IL_T_dB,IL_T_dB_M))
title('IL(\omega)')
xlim([25 6000])

%%
figure(4)
semilogx(freq,10*log10(P_sans2./meas_D),'+-')
hold on
semilogx(Tc_double(:,1),10*log10(1./Tc_double(:,2)),'-')
semilogx(Tc_double(:,1),10*log10(1./Tc_double(:,3)),'-')


xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
xlim([25 18000])
legend(sprintf('Double, IL = %.1f dB, IL_M = %.1f dB',IL_D_dB,IL_D_dB_M),'Diffus simu','Diffus simu + poreux')
title('IL(\omega)')
xlim([25 6000])

figure(5)
semilogx(freq,10*log10(P_sans2./meas_T),'+-')

hold on
semilogx(Tc_double(:,1),10*log10(1./Tc_triple(:,2)),'-')
semilogx(Tc_double(:,1),10*log10(1./Tc_triple(:,3)),'-')
xlabel('Frequency [Hz]')
ylabel('Amplitude [dB]')
xlim([25 18000])
legend(sprintf('Triple, IL = %.1f dB, IL_M = %.1f dB',IL_T_dB,IL_T_dB_M),'Diffus simu','Diffus simu + poreux')
title('IL(\omega)')
xlim([25 6000])



%  FigurePlacecement(1)
