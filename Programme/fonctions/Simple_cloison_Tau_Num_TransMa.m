function [ Tau ] = Simple_cloison_Tau_Num_TransMa(h,E,f,k,theta,rho,nu,Zc1,Zc2,eta)
% [ Tau ] = Simple_cloison_Tau_Num_TransMa(h,E,f,k,theta,rho,nu,Zc1,Zc2,eta)
% h1 épaisseur plaque 1
% E1 module d'young plaque 1
% f frequency (vector (1x1xf))
% k nombre d'onde vecteur d'entrée
% theta
% rho1  densités plaque 1  
% nu
% Zc1 impedance entree
% Zc2 impedance sortie
% eta facteur amortissement plaque
if nargin<10
eta=0;
end
E=E.*(1+1j*eta);
mu= rho*h; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E*h^3/(12*(1-nu^2));

%%
kp = (omega).^(1/4).*sqrt(mu/D);

%%
I = ones(1,1,length(f));
o = zeros(1,1,length(f));

Zpl = D*(k.^4.*sin(theta).^4-kp.^(2).*omega.^(3/2))./( 1j.*omega); 

Tpl = [ I Zpl ; o I];
Tau = 4.*real(Zc1)./real(Zc2) .* abs(Tpl(1,1,:)+Tpl(1,2,:)./Zc2.*cos(theta)+Zc1./cos(theta).*Tpl(2,1,:)+Tpl(2,2,:).*Zc1./Zc2).^(-2);
Tau = permute(Tau,[ 3 2 1]);
end