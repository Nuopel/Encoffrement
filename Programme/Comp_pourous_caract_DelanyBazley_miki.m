% 
%  File    :   Test_Delany-Bazley-Miki.m 
% 
% Author1: Samuel Dupont
% Date:    January 2017 
%   
%  Description :  
%                     
%    1) This program compare Delany-Bazley and Miki coefficients incidence normal.        
%                     
%    2) It use vect
%                     

%
%  Examples of how to utilize the script: The characteristics of the porous
%     material can be change in the first section

clear variables; close all; clc ;
addpath(genpath('../Toolbox'));

%% Define the different parameters
    f = [20:20000].';
    omega = 2 * pi .*f;
   
    rho_0 = 1.213;      % [Kg.m-3] density at rest of air at 18C, 1atm
    c_0   = 342.2;      % [m.s-1] speed of sound in air at 18C, 1atm
    sigma = 30000 ;       % [N.s.m-4] static air flow resistivity of material
    h     = 0.05 ;       % [m] thickness of material
    phi   = 0.98 ;
    lambda = 70 ; 
    lambdap  = 210;
    tortu = 1.1; 

%% Verification
    fprintf('The valid range for the Porous material is %g-%g Hz \n',0.01 * sigma, sigma) 
    U = rho_0 * f / sigma;

%% Delany and Bazley model
    
    Z_DB = rho_0*c_0 * ( 1 + 0.0571*U.^(-0.754) - 1i*0.087*U.^(-0.732) ); 
    km_DB = omega./c_0 .* ( 1 + 0.0978*U.^(-0.7) - 1i * 0.189*U.^(-0.595) );
    
%% Miki Model
    Z_MIK = rho_0*c_0 * ( 1 + 0.0785*U.^(-0.632) - 1i*0.120*U.^(-0.632) ); 
    km_MIK = omega./c_0 .* ( 1 + 0.122*U.^(-0.618) - 1i * 0.180*U.^(-0.618) );
    
    
%%  Calculation for a rigid Wall
    Z_0 = rho_0*c_0;

    % Delany Bazley
    Z_p_DB = -1i.* Z_DB .*cot(km_DB*h);
    alpha_DB = 1 - ( abs( (Z_p_DB-Z_0)./ (Z_p_DB+Z_0) ) ).^2;
    % Miki
    Z_p_MIK = -1i.*Z_MIK .*cot(km_MIK*h);
    alpha_MIK = 1 - ( abs( (Z_p_MIK-Z_0)./( Z_p_MIK +Z_0) ) ).^2;
    
%% Plot
    figure(1)
    subplot(211)
        plot(f, real([Z_p_DB Z_p_MIK])/(rho_0*c_0))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Re(Z)')
        legend('Delany Bazley','Miki')
    
    subplot(212)
        plot(f, imag([Z_p_DB Z_p_MIK])/(rho_0*c_0))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Imag(Z)')
        legend('Delany Bazley','Miki')
    
    figure(2)
        plot(f, ([alpha_DB alpha_MIK]))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Absorption ')
        legend('Delany Bazley','Miki')
