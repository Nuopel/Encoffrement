clear variables; clc; close all
addpath(genpath('./fonctions'));

f=20:10:20e3;
omega = 2*pi*f;
%% parameters of the plate
h1=5e-3; % epaisseur de la plaque 1
h2=10e-3; % epaisseur de la plaque 2

eta = 0.01;
E1=12.6e9*(1+eta*1j); %young modulus plaque 1 
E2=13.6e9*(1+eta*1j); %young modulus plaque 1 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid surfacique de la plaque 1 

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));

e=0.1;
%% parameters of the porous material
    %melamine
%     sigma =10900;      % [N.s.m-4] Static air flow resistivity of material
%     h     = 2*2.54e-2;       % [m] Thickness of material
%     phi_p   = 0.99;  % [/] Porosity
%     lambda = 100e-6 ;     % [um] Viscous length
%     lambdap  = 130e-6;    % [um] Thermic length 
%     tortu = 1;    % [/] Tortuosity
%     
%     %% glass wool ?
    sigma =9000;      % [N.s.m-4] Static air flow resistivity of material
    h     = e;       % [m] Thickness of material
    phi_p   = 0.968;  % [/] Porosity
    lambda = 57e-6 ;     % [um] Viscous length
    lambdap  = 123e-6;    % [um] Thermic length 
    tortu = 1.0295;    % [/] Tortuosity
%     
%         %% rock wool ?
%     sigma =20600;      % [N.s.m-4] Static air flow resistivity of material
%     h     = e;       % [m] Thickness of material
%     phi_p   = 0.98;  % [/] Porosity
%     lambda = 120e-6 ;     % [um] Viscous length
%     lambdap  = 128e-6;    % [um] Thermic length 
%     tortu = 1.01;    % [/] Tortuosity

% [Zp,kp,rhof,Keff] = ChampouxA1j_coef_v2(omega,phi_p,sigma,tortu,lambda,lambdap);
% cf = Zp./rhof;
%% parameters of the air 

rho0=1.2;

thetad=60;
theta=thetad*pi/180;

%% calcul pour cloison simple

Tau_p1 = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule
TauD =    Double_cloison_Tau_Num(h1,h2,E1,E2,permute(f,[1 3 2]),theta,rho,nu,nu,e);

%% calcul amortissement
sigma = 4000;
[ Zf,kp ] = DelanyBazleyMiki_Coef( sigma, f );%% Plot zone

    %  Calculation for a rigid Wall
    Z_0 = 343*1.2;
    % Miki
    Z_p_MIK = -1i.*Zf .*cot(kp*e);
    alpha_MIK = 1 - ( abs( (Z_p_MIK-Z_0)./( Z_p_MIK +Z_0) ) ).^2;
       
    rhof  = Zf.*kp ./ (omega);
    cf = Zf./rhof;

%% calcul plaque + armotissement
Tau_p1_A  = Simple_cloison_Amorti_Tau_Num(h1,E1,f,theta,rho,nu,rhof,cf);
TauD_A = Double_cloison_Tau_amortissement_Num(h1,h2,E1,E2,permute(f,[1 3 2]),theta,rho,nu,nu,e,permute(rhof,[1 3 2]),permute(cf,[1 3 2]),permute(kp,[1 3 2]));


    
%% Plot

c0=343;
fc1=c0^2/2/pi/(sin(theta)^2)*sqrt(mu1/D1);% fréquence de coincidence p1
fc2=c0^2/2/pi/(sin(theta)^2)*sqrt(mu2/D2);% fréquence de coincidence p2fc1

    figure(1)
    subplot(211)
        plot(f, real( Z_p_MIK)/Z_0)
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Re(Z)')
        legend('Miki')
    
    subplot(212)
        plot(f, imag( Z_p_MIK)/Z_0)
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Imag(Z)')
        legend('Miki')
    
    figure(2)
        plot(f, ( alpha_MIK))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Absorption ')
        legend('Miki')
    
    figure(3)
        semilogx(f,10*log10(1./abs(Tau_p1)))
        hold  on
        semilogx(f,10*log10(1./abs(Tau_p1_A)))
        xlabel('Frequence [Hz] log')
        ylabel('Indice d''affaiblissement 10log_{10}')
        legend('Plaque1','Plaque1+amortissement')
        xlim([f(1) f(end)])

    figure(4)
        semilogx(f,10*log10(1./abs(TauD)))
        hold  on
        semilogx(f,10*log10(1./abs(TauD_A)))
        xlabel('Frequence [Hz] log')
        ylabel('Indice d''affaiblissement 10log_{10}')
        for ii = 1:5
        plot([1*c0*ii/(2*e*cosd(thetad)) 1*c0*ii/(2*e*cosd(thetad))],[0 250],'g')
        end
        xlim([f(1) f(end)])
        plot([fc1 fc1],[0 250],'r')
        plot([fc2 fc2],[0 250],'r')
        plot([1*c0/(2*e*cosd(thetad)) 1*c0/(2*e*cosd(thetad))],[0 250],'g')
        xlim([f(1) f(end)])
            legend('cloison','cloison+amortissement')

    FigurePlacecement(1)