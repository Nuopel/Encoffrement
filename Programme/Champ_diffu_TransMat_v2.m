clear variables; clc; close all
addpath(genpath('./fonctions'));

%% parameters of the plate
h1=10e-3; % epaisseur de la plaque 1
h2=20e-3; % epaisseur de la plaque 1
h3 =15e-3;
% 
E1=12.6e9; %young modulus plaque 1 
E2= 12.2e9; %young modulus plaque 2 
E3 = E1;
% 
rho=750; % densité plaque
nu=0.245; % coefficient de poisson MDF
% 
mu1= rho*h1; % poid surfacique de la plaque 1 
mu2= rho*h2; % poid de la plaque 2
mu3= rho*h3; % poid de la plaque 2
% 
D1= E1*h1^3/(12*(1-nu^2));
D2= E2*h2^3/(12*(1-nu^2));
D3= E3*h3^3/(12*(1-nu^2));


e=0.2; % espace entre les cloisons 1 et 2
L=0.1; % espace entre les cloisons 2 et 3

%% porous material
syms theta f
thetalim= 70*pi/180;
c0=343; % vitesse de l'air
rho0=1.2;
omega = 2*pi*f;
k0 = omega/c0;
Z0=rho0*c0;
%   glass wool ?
    sigma =9000;      % [N.s.m-4] Static air flow resistivity of material
    h     = e;       % [m] Thickness of material
    phi_p   = 0.968;  % [/] Porosity
    lambda = 57e-6 ;     % [um] Viscous length
    lambdap  = 123e-6;    % [um] Thermic length 
    tortu = 1.0295;    % [/] Tortuosity

[Zf,kf,rhof,Keff] = ChampouxA1j_coef_v2(omega,phi_p,sigma,tortu,lambda,lambdap);
cf = Zf./rhof;

%% simple cloison

mu= rho*h1; % poid surfacique de la plaque 1 
D= E1*h1^3/(12*(1-nu^2));
eta=0;
Tau =Simple_cloison_Tau_syms_TransMa(h1,E1,f,k0,theta,rho,nu,Z0,Z0,eta); % plaque 1 seule
Taup =Simple_cloison_Tau_syms_AM_TransMa(h1,E1,f,k0,theta,rho,nu,0.02,kf,Zf,Z0,Z0,eta); % plaque 1 seule

%%
fn=20:5:16000;
Tc=zeros(length(fn),2);
T_dif = matlabFunction(Tau*cos(theta)*sin(theta));
T_difp = matlabFunction(Taup*cos(theta)*sin(theta));

tic
%%

mcf=matlabFunction(Tau);
Tc_test=zeros(length(fn),1);
T_diftest = matlabFunction(Tau*cos(theta)*sin(theta));
for N = [100 500 2000 ]
thetaint=linspace(0,70*pi/180,N);
% % 
% Tc=zeros(length(fn),2);
T_diftest = matlabFunction(Tau*cos(theta)*sin(theta));
% %%
for jj= 1:length(thetaint)
    for ii=1:length(fn)
        Tc_test(ii) =Tc_test(ii)+2*feval(T_diftest,fn(ii),thetaint(jj));
    end
    fprintf('%i pc \n',round(jj/length(thetaint)*100))

end
Tc_test =Tc_test./N.*70*pi/180;
semilogx(fn,10*log10(1./Tc_test))
hold on
end
% for ii=1:length(fn)
% Tc(ii,1)=2*integral(@(theta)T_dif(fn(ii),theta),0,thetalim);
% % Tc(ii,2)=2*integral(@(theta)T_difp(fn(ii),theta),0,thetalim);
% fprintf('%i pc \n',round(ii/length(fn)*100))
% end
% toc
% semilogx(fn,10*log10(1./Tc))
% hold on 
% semilogx(fn,10*log10(1./Tc_test))

%% Double cloison
% Tau =Double_cloison_Tau_syms_MatTra(h1,h2,E1,E2,f,k0,theta,rho,rho,nu,nu,e,k0,Z0,Z0,Z0); % double cloison
% Tau = Double_cloison_Tau_symf_theta(h1,h2,E1,E2,f,theta,rho,nu,nu,e);
% Tau2 =Double_cloison_Tau_syms_MatTra(h1,h2,E1,E2,f,k0,theta,rho,rho,nu,nu,e,k0,Z0,Z0,Z0); % double cloison
% 
% % Taup =Double_cloison_Tau_syms_AM_MatTra(h1,h2,E1,E2,f,k0,theta,rho,rho,nu,nu,e,0.02,kf,kf,Zf,Zf,Z0,Z0,eta); % double cloison amortie
% fn=20:10:16000;
% N=10000;
% thetaint=linspace(0,70*pi/180,N);
% % 
% Tc_test=zeros(length(fn),1);
% 
% % Tc=zeros(length(fn),2);
% T_diftest = matlabFunction(Tau*cos(theta)*sin(theta));
% T_difp = matlabFunction(Tau2*cos(theta)*sin(theta));
% % %%
% for jj= 1:length(thetaint)
%     for ii=1:length(fn)
%         Tc_test(ii) =Tc_test(ii)+2*feval(T_diftest,fn(ii),thetaint(jj));
%     end
%     fprintf('%i pc \n',round(jj/length(thetathetaint)*100))
% 
% end
% Tc_test =Tc_test./N.*70*pi/180;
% semilogx(fn,10*log10(1./Tc_test))

% tic
% for ii=1:length(fn)
% Tc(ii,1)=2*integral(@(theta)T_dif(fn(ii),theta),0,thetalim);
% Tc(ii,2)=2*integral(@(theta)T_difp(fn(ii),theta),0,thetalim);
% fprintf('%i pc \n',round(ii/length(fn)*100))
% end
% toc
% semilogx(fn,10*log10(1./abs(Tc)))
% hold on 
%  
% %% Triple cloison
% Tau = triple_cloison_Tau_syms_TransMat(h1,h2,h3,E1,E2,E3,f,k0,theta,rho,rho,rho,nu,nu,nu,e,L,k0,k0,Z0,Z0,Z0,Z0,eta); % triple cloison
% Taup =triple_cloison_Tau_syms_AM_TransMat(h1,h2,h3,E1,E2,E3,f,k0,theta,rho,rho,rho,nu,nu,nu,e,L,0.02,kf,kf,kf,Zf,Zf,Zf,Z0,Z0,eta); % double cloison amortie
% fn=20:10:16000;
% 
% Tc=zeros(length(fn),2);
% T_dif = matlabFunction(Tau*cos(theta)*sin(theta));
% T_difp = matlabFunction(Taup*cos(theta)*sin(theta));
% 
% tic
% for ii=1:length(fn)
% Tc(ii,1)=2*integral(@(theta)T_dif(fn(ii),theta),0,thetalim);
% Tc(ii,2)=2*integral(@(theta)T_difp(fn(ii),theta),0,thetalim);
% fprintf('%i pc \n',round(ii/length(fn)*100))
% end
% toc
% semilogx(fn,10*log10(1./Tc))
% legend('p1','p1 poreux','d','d poreux','t','t poreux')
% xlim([fn(1) fn(end)])
