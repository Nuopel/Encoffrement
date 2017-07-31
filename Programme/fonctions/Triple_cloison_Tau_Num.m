function [ Tau ] = Triple_cloison_Tau_Num(h1,h2,h3,E1,E2,E3,f,thetad,rho,nu1,nu2,nu3,e,L)
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


Tau = zeros(length(f),1);
omega = f*2*pi;
k0 = omega/c0;
k = k0 *sind(thetad);
TH=cosd(thetad);% define constant for matrix
aep = exp(+1j*k0*TH*e);
aem = exp(-1j*k0*TH*e);

ALp = exp(+1j*k0*TH*L);
ALm = exp(-1j*k0*TH*L);

x=1j*k0*TH;
o=zeros(1,1,length(f));
I=ones(1,1,length(f));

A= [  x o o rho0*omega.^2 o o o o o; ...
    o -x x  rho0*omega.^2 o o o o o; ...
    I -I -I -mu1*omega.^2+D1*k.^4 o o o o o;...
    
    o -x.*aem +x.*aep o o o  rho0*omega.^2 o o;...
    o o o o -x.*aem +x.*aep  rho0*omega.^2 o o;...
    o aem +aep  o -aem -aep  -mu2*omega.^2+D2*k.^4 o o;...
    
    o o o o -x.*ALm x.*ALp o o rho0*omega.^2 ;...
    o o o o o o  o -x.*ALm rho0*omega.^2  ;
    o o o o ALm ALp o -ALm -mu3*omega.^2+D3*k.^4
    ];
B = zeros(9,1,length(f));
B(1,1,:)=x;
B(3,1,:)=-1;

for ii=1:length(f)
s =  A(:,:,ii)\B(:,:,ii);
Tau(ii) = s(8)^2;
end

end

