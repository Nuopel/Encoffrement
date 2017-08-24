clear all; clc; close all
addpath(genpath('./fonctions'));

f=permute(21:1:16e3,[ 1 3 2]);


%% parameters of the plate

h1=3/4*2.54e-2; % epaisseur de la plaque 1
h2=1/2*2.54e-2; % epaisseur de la plaque 1
h3=1/4*2.54e-2;
% 
E1=2.3e9; %young modulus plaque 1 
E2=3e9; %young modulus plaque 2
E3=3e9; %young modulus plaque 3 
% 
rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF
% 
mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
mu3= rho*h3; % poid de la plaque 2
% 
D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));
D3= E3*h3^3/(12*(1-nu^2));


e=0.03; % espace entre les cloisons 1 et 2
L=0.07; % espace entre les cloisons 2 et 3

%% parameters of the air 

c0=343; % vitesse de l'air
rho0=1.2;

thetad=50;
theta=thetad*pi/180;
TH = cos(theta);
Tau = zeros(length(f),1);
omega = f*2*pi;
k0 = omega/c0;
k = k0 *sind(thetad);

Tau=Triple_cloison_Tau_Num(h1,h2,h3,E1,E2,E3,f,thetad,rho,nu,nu,nu,e,L);

f=permute(f,[3 2 1]);
omega =  2*pi*f;
k0 = omega/c0;

fc1=c0^2/2/pi/(sind(thetad)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sind(thetad)^2)*sqrt(mu2/D2);% fréquence de coincidence p2
fc3=c0^2/2/pi/(sind(thetad)^2)*sqrt(mu3/D3);% fréquence de coincidence p2

f=permute(f,[ 3 2 1]);
Tau_d1 = Double_cloison_Tau_Num(h1,h2,E1,E2,f,theta,rho,nu,nu,e);
Tau_d2 = Double_cloison_Tau_Num(h2,h3,E2,E3,f,theta,rho,nu,nu,L-e);
Tau_d3 = Double_cloison_Tau_Num(h1,h3,E1,E3,f,theta,rho,nu,nu,L);

f=permute(f,[ 3 2 1]);


Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule
Tau_p2 =  Simple_cloison_Tau_Num(h2,E2,f,theta,rho,nu); % plaque 2 seule
Tau_p3 =  Simple_cloison_Tau_Num(h3,E3,f,theta,rho,nu); % plaque 3 seule

semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./real(Tau_p1)))
semilogx(f,10*log10(1./real(Tau_p2)))
semilogx(f,10*log10(1./real(Tau_p3)))

xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}{1/\tau}')
legend('Cloison triple','Plaque1','Plaque2','Plaque3')
xlim([f(1) f(end)])

figure(2)
semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./abs(Tau_d1)))
semilogx(f,10*log10(1./abs(Tau_d2)))
% semilogx(f,10*log10(1./abs(Tau_d3)))

xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement10log_{10}{1/\tau}')
xlim([f(1) f(end)])
plot([fc1 fc1],[0 200],'r')
plot([fc2 fc2],[0 200],'r')
plot([fc3 fc3],[0 200],'r')
legend('Cloison triple','double h1 h2','double h2 h3')

ylim([0 200])
% fc = wc/2/pi
% fn =1*c0/(2*e*cosd(thetad))
% fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))