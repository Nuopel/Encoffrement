clear all; clc; close all
addpath(genpath('./fonctions'));

f=permute(20:10:20e3,[ 1 3 2]);
omega = 2*pi*f;
I = ones(1,1,length(f));
o = zeros(1,1,length(f));

%% Air parameters
c0 = 343; %celerity of air
k = omega/c0; % wave number in air
rho0=1.2;
Z0 = rho0*c0;
%% parameters of the plate
h1=10e-3; % epaisseur de la plaque 1

E1=12.6e9; %young modulus plaque 1 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
D1= E1*h1^3/(12*(1-nu^2));
eta=0;
wc=c0^2*sqrt(mu1/D1);
%% source inflence
thetad=20;
theta = thetad*pi/180;

%% Transfer matrix linked to the plate
    kp1 = (omega).^(1/4).*sqrt(mu1/D1);

    Zpl2 = 1j.*omega.*mu1.*(1-(omega./wc).^2.*sin(theta)^4.*(1+1j.*eta));
    Zpl1 = D1*(k.^4.*sin(theta).^4-kp1.^(2).*omega.^(3/2))./( 1j.*omega); 
    Tpl = [ I Zpl1 ; o I];
    Tpl2 = [ I Zpl2 ; o I];

Tau  = Simple_cloison_Tau_Num_TransMa(h1,E1,f,k,theta,rho,nu,Z0,Z0,eta);
R3=-10*log10(Tau);
%% transmission loss
R =  20*log10(abs(Tpl(1,1,:)+Tpl(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*Tpl(2,1,:)+Tpl(2,2,:)))-6; 
R = permute(R,[3 1 2]);
R2 =  20*log10(abs(Tpl2(1,1,:)+Tpl2(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*Tpl2(2,1,:)+Tpl2(2,2,:)))-6; 
R2 = permute(R2,[3 1 2]);

figure(1)
semilogx(permute(f,[3 2 1]),R)
hold on
semilogx(permute(f,[3 2 1]),R2,'--')
semilogx(permute(f,[3 2 1]),R3,'--')