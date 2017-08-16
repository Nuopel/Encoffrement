function [ Tau ] = triple_cloison_Tau_Num_AM_MatTra(h1,h2,h3,E1,E2,E3,f,k,theta,rho1,rho2,rho3,nu1,nu2,nu3,e,e2,ef,kf,kf2,kff,Zf,Zf2,Zff,Zc1,Zc2,eta)
% [ Tau ] = Double_cloison_Tau_Num_MatTra(h1,h2,E1,E2,f,theta,rho,nu1,nu2,e,kf,Zf,Zc1,Zc2,eta)
% h1 épaisseur plaque 1
% h2 épaisseur plaque 2
% E1 module d'young plaque 1
% E2 module d'young plaque 1
% f frequency (vector (1x1xf))
% k nombre d'onde vecteur d'entrée
% theta
% rho1 2 densités plaque 1 2 
% nu1
% nu2
% e espace entre les deux plaques
% kf nombre d'onde dans le fluide entre les plaques
% Zf impedance surface fluide entre les plaque
% Zc1 impedance entree
% Zc2 impedance sortie
% eta facteur amortissement plaque
[f1,f2,f3]=size(f);
if f3==1 && f1~=1
   f=permute(f,[3 2 1]); 
end
if f3==1 && f2~=1
   f=permute(f,[1 3 2]); 
end
if nargin<27
eta=0;
end
E1=E1.*(1+1j*eta);
E2=E2.*(1+1j*eta);
E3=E3.*(1+1j*eta);


mu1= rho1*h1; % poid surfacique de la plaque 1 
mu2= rho2*h2; % poid de la plaque 2
mu3= rho3*h3; % poid de la plaque 3

D1= E1*h1^3/(12*(1-nu1^2));
D2= E2*h2^3/(12*(1-nu2^2));
D3= E3*h3^3/(12*(1-nu3^2));

omega = f*2*pi;

I = ones(1,1,length(f));
o = zeros(1,1,length(f));

%% Transfer matrix linked to the porous eq fluid layer in front
thetaf = asin(k.*sin(theta)./kff);
kx = kff.*cos(thetaf);
Tff = [cos(kx*ef) 1j.*Zff./cos(thetaf).*sin(kx.*ef); 1j.*cos(thetaf)./Zff.*sin(kx.*ef) cos(kx*e) ];

%% Transfer matrix linked to the porous eq fluid layer 1
thetaf = asin(k.*sin(theta)./kf);
kx = kf.*cos(thetaf);
Tf = [cos(kx*e) 1j.*Zf./cos(thetaf).*sin(kx.*e); 1j.*cos(thetaf)./Zf.*sin(kx.*e) cos(kx*e) ];

%% Transfer matrix linked to the porous eq fluid layer 2
thetaf2 = asin(k.*sin(theta)./kf2);
kx2 = kf2.*cos(thetaf2);
Tf2 = [cos(kx2*e2) 1j.*Zf2./cos(thetaf2).*sin(kx2.*e2); 1j.*cos(thetaf2)./Zf2.*sin(kx2.*e2) cos(kx2*e2) ];

%% Transfer matrix linked to the plates
    kp1 = (omega).^(1/4).*sqrt(mu1/D1);
    kp2 = (omega).^(1/4).*sqrt(mu2/D2);
    kp3 = (omega).^(1/4).*sqrt(mu3/D3);
       

      Zpl1 =D1*(k.^4.*sin(theta).^4-kp1.^(2).*omega.^(3/2))./( 1j.*omega); 
      Zpl2= D2*(k.^4.*sin(theta).^4-kp2.^(2).*omega.^(3/2))./( 1j.*omega); 
      Zpl3= D3*(k.^4.*sin(theta).^4-kp3.^(2).*omega.^(3/2))./( 1j.*omega); 

    Tpl1 = [ I Zpl1 ; o I];
    Tpl2 = [ I Zpl2 ; o I];
    Tpl3 = [ I Zpl3 ; o I];

%% Global tranfert matrix
T = zeros(2,2,length(f));
for ii = 1:length(f)
    T(:,:,ii)  = Tff(:,:,ii)*Tpl1(:,:,ii)*Tf(:,:,ii)*Tpl2(:,:,ii)*Tf2(:,:,ii)*Tpl3(:,:,ii);
end
%% Tau
Tau = 4.*real(Zc1)./real(Zc2) .* abs(T(1,1,:)+T(1,2,:)./Zc2.*cos(theta)+Zc1./cos(theta).*T(2,1,:)+T(2,2,:).*Zc1./Zc2).^(-2);
Tau = permute(Tau,[ 3 2 1]);

end

