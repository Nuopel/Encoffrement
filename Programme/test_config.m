% clear variables; clc; close all
addpath(genpath('./fonctions'));

f=permute(20:10:20e3,[1 3 2]);
omega = 2*pi*f;
%% parameters of the plate
h1=19e-3; % epaisseur de la plaque 1
h2=12.7e-3; % epaisseur de la plaque 2
h3=6.3e-3; % epaisseur de la plaque 3

E1=2.3e9; %young modulus plaque 1 
E2=3.0e9; %young modulus plaque 2
E3=3.0e9; %young modulus plaque 3 

E22=2.0e9; %young modulus plaque 2

rho=750; % densité plaque
nu=0.25; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid surfacique de la plaque 2 
mu3= rho*h3; % poid surfacique de la plaque 3 

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));
D3= E3*h3^3/(12*(1-nu^2));

e=0.12;
L=0.04;
eta=0.001;
%% parameters of the porous material
%     %% glass wool ?
    sigma =9000;      % [N.s.m-4] Static air flow resistivity of material
    h     = e;       % [m] Thickness of material
    phi_p   = 0.968;  % [/] Porosity
    lambda = 57e-6 ;     % [um] Viscous length
    lambdap  = 123e-6;    % [um] Thermic length 
    tortu = 1.0295;    % [/] Tortuosity

[Zf,kf,rhof,Keff] = ChampouxA1j_coef_v2(omega,phi_p,sigma,tortu,lambda,lambdap);
cf = Zf./rhof;

%% parameters of the air 
rho0=1.2;
c0=343;
Z0=c0*rho0;
k0=omega./c0;
thetad=50;
theta=thetad*pi/180;

%% calcul pour cloison simple
Tau_p1 =  Simple_cloison_Tau_Num_TransMa(h1,E1,f,k0,theta,rho,nu,Z0,Z0,eta); % plaque 1 seule
Tau_p2 =  Simple_cloison_Tau_Num_TransMa(h2,E2,f,k0,theta,rho,nu,Z0,Z0,eta); % plaque 1 seule
Tau_p3 =  Simple_cloison_Tau_Num_TransMa(h3,E3,f,k0,theta,rho,nu,Z0,Z0,eta); % plaque 1 seule

%% calcul plaque double + armotissement
% Tau_p1_A  = Simple_cloison_Tau_Num_TransMa(h1,E1,f,kf,theta,rho,nu,Zf,Z0,eta); % plaque 1 seule
Tau_p1_A  = Simple_cloison_Tau_Num_AM_TransMa(h1,E1,f,k0,theta,rho,nu,Z0,Z0,e,kf,Zf,eta); % plaque 1 seule

Tau_p2_A  = Simple_cloison_Tau_Num_TransMa(h2,E2,f,kf,theta,rho,nu,Zf,Z0,eta); % plaque 1 seule
Tau_p3_A  = Simple_cloison_Tau_Num_TransMa(h3,E3,f,kf,theta,rho,nu,Zf,Z0,eta); % plaque 1 seule

TauD_A =  Double_cloison_Tau_Num_AM_MatTra(h1,h2,E1,E2,f,k0,theta,rho,rho,nu,nu,e,0.02,kf,Zf,Z0,Z0,eta);
TauD=  Double_cloison_Tau_Num_MatTra(h1,h2,E1,E2,f,k0,theta,rho,rho,nu,nu,e,k0,Z0,Z0,Z0,eta);

% TauD2_A =  Double_cloison_Tau_Num_MatTra(h2,h3,E2,E3,f,k0,theta,rho,rho,nu,nu,L,kf,Zf,Z0,Z0,eta);
TauD2_A =  Double_cloison_Tau_Num_AM_MatTra(h2,h3,E2,E3,f,k0,theta,rho,rho,nu,nu,L,0.02,kf,Zf,Z0,Z0,eta);

TauD2 =  Double_cloison_Tau_Num_MatTra(h2,h3,E2,E3,f,k0,theta,rho,rho,nu,nu,L,k0,Z0,Z0,Z0,eta);


%% calcul plaque triple + armotissement
TauT = triple_cloison_Tau_Num_MatTra(h1,h2,h3,E1,E2,E3,f,k0,theta,rho,rho,rho,nu,nu,nu,e,L,k0,k0,Z0,Z0,Z0,Z0,eta);
% TauT_A = triple_cloison_Tau_Num_MatTra(h1,h2,h3,E1,E2,E3,f,k0,theta,rho,rho,rho,nu,nu,nu,e,L,kf,kf,Zf,Zf,Z0,Z0,eta);
TauT_A = triple_cloison_Tau_Num_AM_MatTra(h1,h2,h3,E1,E2,E3,f,k0,theta,rho,rho,rho,nu,nu,nu,e,L,0.03,kf,kf,kf,Zf,Zf,Zf,Z0,Z0,eta);

%% Plot

c0=343;
fc1=c0^2/2/pi/(sin(theta)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sin(theta)^2)*sqrt(mu2/D2);% fréquence de coincidence p2fc1
fc3=c0^2/2/pi/(sin(theta)^2)*sqrt(mu3/D3);% fréquence de coincidence p2fc1
        figure(1)
        semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauT_A)))
        hold on
        figure(2)
        semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauD_A)))
        hold on 
     
        figure(3)
%         semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauT)))
%         hold  on
        semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauT_A)));
        hold on
%         semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauD_A)))
        semilogx(permute(f,[3,2,1]),10*log10(1./abs(TauD_A)))

        xlabel('Frequence [Hz] log')
        ylabel('Indice d''affaiblissement 10log_{10}')
        xlim([f(1) f(end)])
%         plot([fc1 fc1],[0 250],'r')
%         plot([fc2 fc2],[0 250],'r')
%         plot([fc3 fc3],[0 250],'r')

        xlim([f(1) f(end)])
        legend('triple','triple+amortissement','double+amortissement','double+amortissement')

    FigurePlacecement(1)