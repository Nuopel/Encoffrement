clear all; clc; close all
addpath(genpath('./fonctions'));

f=permute(20:1:40e3,[ 1 3 2]);


%% parameters of the plate
h1=35e-3; % epaisseur de la plaque 1
h2=20e-3; % epaisseur de la plaque 1
h3 =10e-3;

E1=12.6e9; %young modulus plaque 1 
E2= 12.2e9; %young modulus plaque 2 
E3 = E1;

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
mu3= rho*h3; % poid de la plaque 2

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));
D3= E3*h3^3/(12*(1-nu^2));

%% parameters of the air 

e=0.1; % espace entre les cloisons 1 et 2
L=0.25; % espace entre les cloisons 2 et 3
c0=343; % vitesse de l'air
rho0=1.2;

thetad=40;
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

Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,thetad,rho,nu); % plaque 1 seule
Tau_p2 =  Simple_cloison_Tau_Num(h2,E2,f,thetad,rho,nu); % plaque 2 seule
Tau_p3 =  Simple_cloison_Tau_Num(h3,E3,f,thetad,rho,nu); % plaque 3 seule

semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./real(Tau_p1)))
semilogx(f,10*log10(1./real(Tau_p2)))
semilogx(f,10*log10(1./real(Tau_p3)))

xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}')
legend('Cloison double','Plaque1','Plaque2','Plaque3')
xlim([f(1) f(end)])
% fc = wc/2/pi
% fn =1*c0/(2*e*cosd(thetad))
% fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))