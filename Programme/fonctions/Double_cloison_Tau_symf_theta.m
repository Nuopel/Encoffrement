function [ Tau ] = Double_cloison_Tau_symf_theta(h1,h2,E1,E2,f,theta,rho,nu1,nu2,e)
syms f theta
c0=343; % vitesse de l'air
rho0=1.2;


mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
D1= E1*h1^3/(12*(1-nu1^2));
D2= E2*h2^3/(12*(1-nu2^2));


omega = f*2*pi;
k0 = omega/c0;
k = k0 *sin(theta);
TH=cos(theta);
alp = exp(+1j*k0*TH*e);
alm = exp(-1j*k0*TH*e);

A= [ 1j*k0*TH 0 0 0 rho0*omega.^2 0; ...
    0 -1j*k0*TH 1j*k0*TH  0 rho0*omega.^2 0; ...
    1 -1 -1 0 -mu1*omega.^2+D1*k.^4 0 ;...
    0 -1j*k0*TH.*alm 1j*k0*TH.*alp 0 0  rho0*omega.^2;...
    0 0 0 -1j*k0*TH.*alm 0 rho0.*omega.^2;...
    0 alm alp -alm 0 -mu2*omega.^2+D2*k.^4 ...
    ];
B = [1j*k0*TH ; 0 ;-1; 0; 0 ; 0];
s =  A\B;
Tau = s(4)^2;

end

