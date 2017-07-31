function [ Tau ] = Simple_cloison_Tau_Num(h,E,f,theta,rho,nu)

c0=343; % vitesse de l'air
rho0=1.2;
TH=cos(theta);
mu= rho*h; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E*h^3/(12*(1-nu^2));
k0 = omega/c0;
Tau = omega.^2 *(rho0*c0)^2/TH^2*4./abs(-omega.^2*mu+D*k0.^4*sin(theta)^4+2*rho0*omega*c0/1j/TH).^2; % plaque 1 seule


end

