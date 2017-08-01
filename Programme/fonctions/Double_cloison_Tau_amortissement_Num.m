function [ Tau ] = Double_cloison_Tau_amortissement_Num(h1,h2,E1,E2,f,theta,rho,nu1,nu2,e,rhop,cp,kp)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

rho0=1.2;
c0=343;

mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
D1= E1*h1^3/(12*(1-nu1^2));
D2= E2*h2^3/(12*(1-nu2^2));


Tau = zeros(length(f),1);
omega = f*2*pi;
k0 = omega./c0;
k = k0 *sin(theta);

thetap = asin(k0.*sin(theta)./kp);
THp=cos(thetap);
alpp = exp(+1j*kp.*THp*e);
almp = exp(-1j*kp.*THp*e);


TH=cos(theta);
alp = exp(+1j*k0.*TH*e);
alm = exp(-1j*k0.*TH*e);

o=zeros(1,1,length(f));
I=ones(1,1,length(f));
A= [ 1j*k0*TH o o o rho0.*omega.^2 o; ...
    o -1j*kp.*THp 1j*kp.*THp o rhop.*omega.^2 o; ...
    I -I -I o -mu1*omega.^2+D1*k.^4 o ;...
    o -1j*kp.*THp.*almp 1j*kp.*THp.*alpp o o  rhop.*omega.^2;...
    o o o -1j*k0.*TH.*alm o rho0.*omega.^2;...
    o almp alpp -alm o -mu2*omega.^2+D2*k.^4 ...
    ];
B = [1j*k0*TH ; o ;-I; o; o;o];
for ii=1:length(f)
s =  A(:,:,ii)\B(:,:,ii);
Tau(ii) = s(4)^2;
end

end

