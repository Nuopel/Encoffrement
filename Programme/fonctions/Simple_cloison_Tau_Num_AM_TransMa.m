function [ Tau ] = Simple_cloison_Tau_Num_AM_TransMa(h,E,f,k,theta,rho,nu,Zc1,Zc2,e,kf,Zf,eta)
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
if nargin<13
eta=0;
end
E=E.*(1+1j*eta);
mu= rho*h; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E*h^3/(12*(1-nu^2));

%%
kp = (omega).^(1/4).*sqrt(mu/D);
%% Transfer matrix linked to the porous eq fluid layer
thetaf = asin(k.*sin(theta)./kf);
kx = kf.*cos(thetaf);
Tf = [cos(kx*e) 1j.*Zf./cos(thetaf).*sin(kx.*e); 1j.*cos(thetaf)./Zf.*sin(kx.*e) cos(kx*e) ];

%%
I = ones(1,1,length(f));
o = zeros(1,1,length(f));

Zpl = D*(k.^4.*sin(theta).^4-kp.^(2).*omega.^(3/2))./( 1j.*omega); 

Tpl = [ I Zpl ; o I];
%% Global tranfert matrix
T = zeros(2,2,length(f));
for ii = 1:length(f)
    T(:,:,ii)  =  Tf(:,:,ii)*Tpl(:,:,ii);
end
%% Tau
Tau = 4.*real(Zc1)./real(Zc2) .* abs(T(1,1,:)+T(1,2,:)./Zc2.*cos(theta)+Zc1./cos(theta).*T(2,1,:)+T(2,2,:).*Zc1./Zc2).^(-2);
Tau = permute(Tau,[ 3 2 1]);
end