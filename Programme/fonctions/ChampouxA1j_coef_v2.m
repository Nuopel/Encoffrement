function [Zcf,kf,Rhoeff,Keff] = ChampouxA1j_coef_v2(om,phi,sig,tortu,lambda,lambdap)
% [Zcf,kf] = ChampouxA1j_coef(om,phi,sig,tortu,lamb,lambp,k0p)

% Author1: Samuel Dupont
% Date:    1january 2017
% 
% Function   : ChampouxA1j_coef 
% 
% Description: Generate Champoux Allard 1johnson coefficient  for pourous
%              material acoording to the book "Propagation of sound in porous 
%              media: modelling sound absorbing materials", page 90
% 
% Parameters : input          
%                             omega,    Q angular frequency
%                                       example = 2 * pi * (20:2000);
% 
%                             phi,      Porosity of the medium
%                                       example = 0.95 ;
% 
%                             sig,      [N.s.m-4] static air flow resistivity of material
%                                       example = 10000 ; 
% 
%                             tortu,    une tortue
%                                       example  = franklin la tortue
% 
%                             lamb,     Length Viscous coefficient
%                                       example = 0.0002;
% 
%                             lamb',    Length Thermic coefficient
%                                       example = 0.0005;
% 
% 
% Return     :    Zcf    Impedance caracteristique du milieu
%                 kf     Wave number in the medium
%                
%  
% Examples of Usage: 
% 

%Air parameters
gama=1.4;
r0=1.204;        % DENSITY            (kg.m-3)
eta=  0.184E-04; % VISCOSITY OF FLUID (kg/m/s)
Pr=0.71;          % PRANDTL NUMBER  
P0= 0.10132e6 ;  %  ATM. PRESSURE      (Pa)     


%%
G = sqrt( 1 + 4 * 1i *tortu^2 * eta * r0 * om  ./ (sig^2 * lambda^2 * phi^2)  );
Rhoeff = r0 * tortu * (1 + (sig *phi)./( 1i * tortu .*r0 *om) .*G );

Gp = sqrt( 1 + ( 1i * r0 * lambdap^2 * Pr * om ) ./ (16 * eta)   );
Keff = gama * P0 ./ ( gama - ( gama - 1 ).*( 1 + (8 * eta *Gp)./(1i*lambdap^2*r0*Pr*om)).^-1);

kf=om.*sqrt(Rhoeff./Keff);
Zcf=sqrt(Rhoeff.*Keff);

end
