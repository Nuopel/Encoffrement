clear all; clc; close all


addpath(genpath('./fonctions'));

f=20:5:16e3;


%% parameters of the plate
h1=1/2*2.54e-2; % epaisseur de la plaque 1
E1=3.6e9; %young modulus plaque 1 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 

D1= E1*h1^3/(12*(1-nu^2));
c0=343;

%% parameters of the air 
thetad=40;
theta = thetad*pi/180;

fc1=c0^2/2/pi/(sin(theta)^2)*sqrt(mu1/D1);% fréquence de coincidence p1

Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule


%% calculate for a finite sized plate
k0=2*pi*f./c0;

a=0.4;% size of the pannel height
b=a;% size of the pannel width
L=sqrt(a*b);

kp = (2*pi*f).^(1/2).*(mu1/D1).^(1/4);
F2integralfn = @(kr,kpi,k0i) sin((kr-kpi).*L/2).^2./( ((kr-kpi).*L/2).^2.*sqrt(k0i.^2-kr.^2) );
Tc=zeros(length(f),1);
for ii = 1:length(f)
Tc(ii,1)=integral(@(kr)F2integralfn(kr,kp(ii),k0(ii)),0,k0(ii));
fprintf('%i pc \n',round(ii/length(f)*100))
end
sig = L.*k0./2/pi.*Tc.';
sig_inf = 1./cos(theta);
Tau_finite = Tau_p1.*(sig./sig_inf).^2;
%% Plot

semilogx(f,10*log10(1./real(Tau_p1)))
hold on
semilogx(f,10*log10(1./real(Tau_finite)))
xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}1/\tau')
plot([fc1 fc1],[0 250],'r')
xlim([f(1) f(end)])
legend('Cloison simple','cloison simple finie')
ylim([0 100])
figure 
semilogx(f,20*log10(abs(sig)))
hold on
semilogx(f,20*log10(abs(sig_inf)))
xlim([f(1) f(end)])

xlabel('Frequence [Hz] log')
ylabel('Efficacité de rayonnement \sigma dB')

% % fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))