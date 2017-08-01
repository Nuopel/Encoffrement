function [ Z,km ] = DelanyBazleyMiki_Coef( sigma, f, opt )
% [ Z_MIK,km_MIK ] = DelanyBazleyMiki_Coef(sigma, h,f )
% 
% Author1: Samuel Dupont
% Date:    January 2017
% 
% Function   : DelanyBazleyMiki_Coef 
% 
% Description: Generate Delany-Bazley coefficient modified by Miki
% 
% Parameters : input          sigma, example = 10000 ;       % [N.s.m-4] static air flow resistivity of material
%                             f,     example = 20:2000
%                             opt,  option if == 1 choose Delany-Bazley
%                                   instead of Miki

% Return     :    Z_MIK Impedance caracteristique du milieu
%                 km_MIK wav number in the medium
% 
% Examples of Usage: 
% 

    if nargin < 3
        opt = 0;
    end

%% define constant
    rho_0 = 1.213;      % [Kg.m-3] density at rest of air at 18C, 1atm
    c_0   = 342.2;      % [m.s-1] speed of sound in air at 18C, 1atm
    omega = 2 * pi * f; % [s-1] angular frequencies

%% Verification
    fprintf('The valid range for the Porous material is %g-%g Hz \n',0.01 * sigma, sigma) 
    U = rho_0 * f / sigma;
    
%% Coefficients
    if opt == 0
        Z = rho_0*c_0 * ( 1 + 0.0785*U.^(-0.632) - 1i*0.120*U.^(-0.632) ); 
        km = omega./c_0 .* ( 1 + 0.122*U.^(-0.618) - 1i * 0.180*U.^(-0.618) );
    else
        Z = rho_0*c_0 * ( 1 + 0.0571*U.^(-0.754) - 1i*0.087*U.^(-0.732) ); 
        km = omega./c_0 .* ( 1 + 0.0978*U.^(-0.7) - 1i * 0.189*U.^(-0.595) );
    end
end