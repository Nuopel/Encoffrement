clear all; clc; close all
addpath(genpath('./fonctions'));

f=permute(20:1:30e3,[ 1 3 2]);
omega = 2*pi*f;
%% Air parameters
c0 = 343; %celerity of air
k = omega/c0; % wave number in air
rho0=1.2;
Z0=c0*rho0;
I = ones(1,1,length(f));
o = zeros(1,1,length(f));

%% parameters of the plate
h1=5e-3; % epaisseur de la plaque 1
h2=10e-3; % epaisseur de la plaque 2
eta=0.0;
E1=12.6e9.*(1+1j*eta); %young modulus plaque 1 
E2= 12.2e9.*(1+1j*eta); %young modulus plaque 2 

rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid surfacique de la plaque 1 

D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));

e=0.05; % space in etween the plates

%% source inflence
thetad=20;
theta = thetad*pi/180;



%% Transfer matrix linked to the air fluid layer
    kx = k*cos(theta);
    Tf = [cos(kx*e) 1j.*omega.*rho0./kx.*sin(kx.*e); 1j.*kx./(omega.*rho0).*sin(kx.*e) cos(kx*e) ];

%% Transfer matrix linked to the porous eq fluid layer
% sigma = 4000;
% [ Zf,kp ] = DelanyBazleyMiki_Coef( sigma, f );%% Plot zone
    % rock wool ?
    sigma =1600;      % [N.s.m-4] Static air flow resistivity of material
    h     = e;       % [m] Thickness of material
    phi_p   = 0.98;  % [/] Porosity
    lambda = 120e-6 ;     % [um] Viscous length
    lambdap  = 128e-6;    % [um] Thermic length 
    tortu = 1.01;    % [/] Tortuosity

[Zf,kp,rhof,Keff] = ChampouxA1j_coef_v2(omega,phi_p,sigma,tortu,lambda,lambdap);
kx = kp*cos(theta);
thetap = asin(k.*sin(theta)./kp);

Tfp = [cos(kx*e) 1j.*Zf./cos(theta).*sin(kx.*e); 1j.*cos(theta)./Zf.*sin(kx.*e) cos(kx*e) ];
%% Transfer matrix linked to the plates
    kp1 = (omega).^(1/4).*sqrt(mu1/D1);
    kp2 = (omega).^(1/4).*sqrt(mu2/D2);

      Zpl1 = D1*(k.^4.*sin(theta).^4-kp1.^(2).*omega.^(3/2))./( 1j.*omega); 
      Zpl2= D2*(k.^4.*sin(theta).^4-kp2.^(2).*omega.^(3/2))./( 1j.*omega); 
    
    Tpl1 = [ I Zpl1 ; o I];
    Tpl2 = [ I Zpl2 ; o I];
%% Global tranfert matrix
T = zeros(2,2,length(f));
Tp = T;

for ii = 1:length(f)
    T(:,:,ii)  =  Tpl1(:,:,ii)*Tf(:,:,ii)*Tpl2(:,:,ii);
    Tp(:,:,ii) =  Tpl1(:,:,ii)*Tfp(:,:,ii)*Tpl2(:,:,ii);

end
%% transmission loss
R =  20*log10(abs(T(1,1,:)+T(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*T(2,1,:)+T(2,2,:)))-6; 
R = permute(R,[3 1 2]);
Rpo =  20*log10(abs(Tp(1,1,:)+Tp(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*Tp(2,1,:)+Tp(2,2,:)))-6; 
Rpo = permute(Rpo,[3 1 2]);


R_p1 =  20*log10(abs(Tpl1(1,1,:)+Tpl1(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*Tpl1(2,1,:)+Tpl1(2,2,:)))-6; 
R_p1 = permute(R_p1,[3 1 2]);

R_p2 =  20*log10(abs(Tpl2(1,1,:)+Tpl2(1,2,:)/Z0*cos(theta)+Z0/cos(theta)*Tpl2(2,1,:)+Tpl2(2,2,:)))-6; 
R_p2 = permute(R_p2,[3 1 2]);

figure(4)
hold on
semilogx(permute(f,[3 2 1]),R)
hold on 
semilogx(permute(f,[3 2 1]),R_p1)
semilogx(permute(f,[3 2 1]),R_p2)
semilogx(permute(f,[3 2 1]),Rpo,'g')

xlim([f(1) f(end)])
legend('Cloison double','Plaque 1', 'Plaque2')


