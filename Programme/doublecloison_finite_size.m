clear all; clc; close all


addpath(genpath('./fonctions'));

f=permute(20:5:16e3,[ 1 3 2]);


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

e=0.05; % espace entre les cloisons
thetad=40;
theta = thetad*pi/180;
[ Tau ] = Double_cloison_Tau_Num(h1,h2,E1,E2,f,theta,rho,nu,nu,e);

%% calculate for a finite sized plate
c0=343;
f=permute(f,[ 3 1 2]);

k0=2*pi*f./c0;
a=0.45;% size of the pannel height
b=0.45;% size of the pannel width
L=sqrt(a*b);

kp = (2*pi*f).^(1/2).*(mu1/D1).^(1/4);
F2integralfn = @(kr,kpi,k0i) sin((kr-kpi).*L/2).^2./( ((kr-kpi).*L/2).^2.*sqrt(k0i.^2-kr.^2) );
Tc=zeros(length(f),1);
for ii = 1:length(f)
Tc(ii,1)=integral(@(kr)F2integralfn(kr,kp(ii),k0(ii)),0,k0(ii));
fprintf('%i pc \n',round(ii/length(f)*100))
end
sig = L.*k0./2/pi.*Tc;
sig_inf = 1./cos(theta);
Tau_finite = Tau.*(sig./sig_inf).^2;
%% Plot
 f=permute(f,[ 3 1 2]);

fc1=c0^2/2/pi/(sin(theta)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sin(theta)^2)*sqrt(mu2/D2);% fréquence de coincidence p2fc1

Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule
Tau_p2 = Simple_cloison_Tau_Num(h2,E1,f,theta,rho,nu);%plaque 2 seule

semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./abs(Tau_finite)))
semilogx(f,10*log10(1./real(Tau_p1)))
semilogx(f,10*log10(1./real(Tau_p2)))
xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}')
plot([fc1 fc1],[0 250],'r')
plot([fc2 fc2],[0 250],'r')
for ii = 1:5
plot([1*c0*ii/(2*e*cosd(thetad)) 1*c0*ii/(2*e*cosd(thetad))],[0 250],'g')
end
% plot([402 402],[0 250],'b')
xlim([f(1) f(end)])
legend('Cloison double','cloison double finitesize','Plaque 1', 'Plaque2')

% % fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))