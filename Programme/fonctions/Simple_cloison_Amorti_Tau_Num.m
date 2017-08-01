function [ Tau ] = Simple_cloison_Amorti_Tau_Num(h,E,f,theta,rho,nu,rhof,cf)

TH=cos(theta);
mu= rho*h; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E*h^3/(12*(1-nu^2));
kf = omega/cf;
Tau = omega.^2 .*(rhof.*cf).^2/TH^2*4./abs(-omega.^2*mu+D*kf.^4*sin(theta)^4+2*rhof.*omega.*cf/1j/TH).^2; % plaque 1 seule

end

