clear all; clc; close all

f=0:0.51:24e3;


%% parameters of the plate
h1=5e-3; % epaisseur de la plaque 1
h2=20e-3; % epaisseur de la plaque 1

E1=12.6e9; %young modulus plaque 1 
E2= 12.2e9; %young modulus plaque 2 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));

%% parameters of the air 

e=0.15; % espace entre les cloisons
c0=343; % vitesse de l'air
rho0=1.2;

thetad=20;
theta=thetad*pi/180;

Tau = zeros(length(f),1);

for ii=1:length(f)
omega = f(ii)*2*pi;
k0 = omega/c0;
k = k0 *sind(thetad);
TH=cosd(thetad);
alp = exp(+1j*k0*TH*e);
alm = exp(-1j*k0*TH*e);
beta = k0*cos(theta) ;
alpha_1 = exp( 1i * k0 * cos(theta) * e) ;
alpha_2 = exp( -1i * k0 * cos(theta) * e);

A= [ 1j*k0*TH 0 0 0 rho0*omega^2 0; ...
    0 -1j*k0*TH 1j*k0*TH 0 rho0*omega^2 0; ...
    1 -1 -1 0 -mu1*omega^2+D1*k^4 0 ;...
    0 -1j*k0*TH*alm 1j*k0*TH*alp 0 0  rho0*omega^2;...
    0 0 0 -1j*k0*TH*alm 0 rho0*omega^2;...
    0 alm alp -alm 0 -mu2*omega^2+D2*k^4 ...
    ];
B = [1j*k0*TH ; 0 ;-1; 0; 0;0];

s =  A\B;
Tau(ii) = s(4)^2;
end

omega =  2*pi*f;
k0 = omega/c0;

fc1=c0^2/2/pi/(sind(thetad)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sind(thetad)^2)*sqrt(mu2/D2);% fréquence de coincidence p2

Tau_p1 = omega.^2 *(rho0*c0)^2/TH^2*4./abs(-omega.^2*mu1+D1*k0.^4*sind(thetad)^4+2*rho0*omega*c0/1j/TH).^2; % plaque 1 seule
Tau_p2 = omega.^2 *(rho0*c0)^2/TH^2*4./abs(-omega.^2*mu2+D2*k0.^4*sind(thetad)^4+2*rho0*omega*c0/1j/TH).^2;%plaque 2 seule

semilogx(f,10*log10(1./real(abs(Tau))))
hold on
semilogx(f,10*log10(1./real(Tau_p1)))
semilogx(f,10*log10(1./real(Tau_p2)))
xlabel('Frequence [Hz] log')
ylabel('Indice d''affaiblissement 10log_{10}')
legend('Cloison double','Plaque')
xlim([f(1) f(end)])
% fc = wc/2/pi
% fn =1*c0/(2*e*cosd(thetad))
% fre=1/(2*pi)*sqrt(2*rho0*(c0)^2/(e*mu*cosd(thetad)^2))