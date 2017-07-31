function [ Tau ] = Triple_cloison_Tau_symf_theta(h1,h2,h3,E1,E2,E3,rho,nu1,nu2,nu3,e,L)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

c0=343; % vitesse de l'air
rho0=1.2;



mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
mu3= rho*h3; % poid de la plaque 2

D1= E1*h1^3/(12*(1-nu1^2));
D2= E2*h2^3/(12*(1-nu2^2));
D3= E3*h3^3/(12*(1-nu3^2));
syms f theta



omega = f*2*pi;
k0 = omega/c0;
k = k0 *sin(theta);
TH=cos(theta);

aep = exp(+1j*k0*TH*e);
aem = exp(-1j*k0*TH*e);

ALp = exp(+1j*k0*TH*L);
ALm = exp(-1j*k0*TH*L);

x=1j*k0*TH;

A= [  x 0 0 rho0*omega.^2 0 0 0 0 0; ...
    0 -x x  rho0*omega.^2 0 0 0 0 0; ...
    1 -1 -1 -mu1*omega.^2+D1*k.^4 0 0 0 0 0;...
    
    0 -x.*aem +x.*aep 0 0 0  rho0*omega.^2 0 0;...
    0 0 0 0 -x.*aem +x.*aep  rho0*omega.^2 0 0;...
    0 aem +aep  0 -aem -aep  -mu2*omega.^2+D2*k.^4 0 0;...
    
    0 0 0 0 -x.*ALm x.*ALp 0 0 rho0*omega.^2 ;...
    0 0 0 0 0 0  0 -x.*ALm rho0*omega.^2  ;
    0 0 0 0 ALm ALp 0 -ALm -mu3*omega.^2+D3*k.^4
    ];
B = [x; 0; -1; 0; 0; 0; 0; 0; 0];

s =  A\B;
Tau = s(8)^2;


end

