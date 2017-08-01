clear all; clc; close all


addpath(genpath('./fonctions'));

f=permute(20:1:40e3,[ 1 3 2]);


%% parameters of the plate
h1=5e-3; % epaisseur de la plaque 1
h2=10e-3; % epaisseur de la plaque 1

E1=12.6e9; %young modulus plaque 1 
E2= 12.2e9; %young modulus plaque 2 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));

%% parameters of the air 

e=0.01; % espace entre les cloisons
thetad=20;
theta = thetad*pi/180;
[ Tau ] = Double_cloison_Tau_Num(h1,h2,E1,E2,f,theta,rho,nu,nu,e);


%% Plot
f= permute(f,[3 2 1]);
% 
c0=343;
fc1=c0^2/2/pi/(sin(theta)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sin(theta)^2)*sqrt(mu2/D2);% fréquence de coincidence p2fc1

Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule
Tau_p2 = Simple_cloison_Tau_Num(h2,E1,f,theta,rho,nu);%plaque 2 seule

semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./real(Tau_p1)))
semilogx(f,10*log10(1./real(Tau_p2)))
xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}')
plot([fc1 fc1],[0 250],'r')
plot([fc2 fc2],[0 250],'r')
for ii = 1:5
plot([1*c0*ii/(2*e*cosd(thetad)) 1*c0*ii/(2*e*cosd(thetad))],[0 250],'g')
end
plot([402 402],[0 250],'b')
xlim([f(1) f(end)])
legend('Cloison double','Plaque 1', 'Plaque2')

% % fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))