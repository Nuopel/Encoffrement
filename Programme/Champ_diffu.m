clear variable; clc; close all
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
L=0.35; % espace entre les cloisons 2 et 3

%% simple cloison

c0=343; % vitesse de l'air
rho0=1.2;
syms theta f% h1 E1 nu rho
% f=50;
mu= rho*h1; % poid surfacique de la plaque 1 
omega = 2*pi*f;
D= E1*h1^3/(12*(1-nu^2));
k0 = omega/c0;
Tau = Simple_cloison_Tau_Num(h1,E1,f,theta,rho,nu); % plaque 1 seule
% T_dif= 2*vpaintegral(Tau*cos(theta)*sin(theta),theta,0,70*pi/180);
%%
fn=20:5:16000;
Tc=zeros(length(fn),1);
T_dif = matlabFunction(Tau*cos(theta)*sin(theta));
tic
for ii=1:length(fn)
 Tc(ii)=integral(@(theta)T_dif(fn(ii),theta),0,pi/2);
% Tc(ii)=double(subs(T_dif,f,fn(ii)));
fprintf('%i pc \n',round(ii/length(fn)*100))
end
toc
semilogx(fn,10*log10(1./Tc))
hold on 
%% Double cloison
Tau = Double_cloison_Tau_symf_theta(h1,h2,E1,E2,f,theta,rho,nu,nu,e);
% T_dif= 2*vpaintegral(Tau*cos(theta)*sin(theta),theta,0,70*pi/180);
%
fn=20:25:16000;
Tc=zeros(length(fn),1);

T_dif = matlabFunction(Tau*cos(theta)*sin(theta));
tic
for ii=1:length(fn)
    Tc(ii)=2*integral(@(theta)T_dif(fn(ii),theta),0,50*pi/180);
    % Tc(ii)=double(subs(T_dif,f,fn(ii)));
    fprintf('%i pc \n',round(ii/length(fn)*100))
end
toc
semilogx(fn,10*log10(abs((1./Tc))))


%% Triple cloison
Taut  = Triple_cloison_Tau_symf_theta(h1,h2,h3,E1,E2,E3,rho,nu,nu,nu,e,L);
% T_dif= 2*vpaintegral(Tau*cos(theta)*sin(theta),theta,0,70*pi/180);
fn=20:25:16000;
T_dift = matlabFunction(Taut*cos(theta)*sin(theta));
Tc=zeros(length(f),1);

tic
for ii=1:length(fn)
Tc(ii)=2*integral(@(theta)T_dift(fn(ii),theta),0,50*pi/180);
% Tc(ii)=double(subs(T_dif,f,fn(ii)));
fprintf('%i pc \n',round(ii/length(fn)*100))
end
toc
semilogx(fn,10*log10(abs(1./Tc)))