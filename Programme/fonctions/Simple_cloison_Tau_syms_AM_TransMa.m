function [ Tau ] = Simple_cloison_Tau_syms_AM_TransMa(h,E,f,k,theta,rho,nu,ef,kf,Zf,Zc1,Zc2,eta)
% [ Tau ] = Simple_cloison_Tau_Num_TransMa(h,E,f,k,theta,rho,nu,Zc1,Zc2,eta)
% h1 épaisseur plaque 1
% E1 module d'young plaque 1
% f frequency (syms)
% k nombre d'onde vecteur d'entrée (syms)
% theta
% rho1  densités plaque 1  
% nu
% Zc1 impedance entree
% Zc2 impedance sortie
% eta facteur amortissement plaque
if nargin<13
eta=0;
end
E=E.*(1+1j*eta);
mu= rho*h; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E*h^3/(12*(1-nu^2));

%% Transfer matrix linked to the porous eq fluid layer in front of the source
thetaf = asin(k.*sin(theta)./kf);
kx = kf.*cos(thetaf);
Tf = [cos(kx*ef) 1j.*Zf./cos(thetaf).*sin(kx.*ef); 1j.*cos(thetaf)./Zf.*sin(kx.*ef) cos(kx*ef) ];
%%
kp = (omega).^(1/4).*sqrt(mu/D);

%%
Zpl = D*(k.^4.*sin(theta).^4-kp.^(2).*omega.^(3/2))./( 1j.*omega); 

Tpl = [ 1 Zpl ; 0 1];
Tpl =Tf*Tpl;
Tau = 4.*real(Zc1)./real(Zc2) .* abs(Tpl(1,1)+Tpl(1,2)./Zc2.*cos(theta)+Zc1./cos(theta).*Tpl(2,1)+Tpl(2,2).*Zc1./Zc2).^(-2);
end