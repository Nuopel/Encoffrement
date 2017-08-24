% 
%  File    :   Comp_Champoux-Allard-Jonh.m 
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

    %         %% rock wool ?
    sigma =20600;      % [N.s.m-4] Static air flow resistivity of material
    h     = 0.1;       % [m] Thickness of material
    phi   = 0.98;  % [/] Porosity
    lambda = 85e-6 ;     % [um] Viscous length
    lambdap  = 90e-6;    % [um] Thermic length 
    tortu = 1.01;    % [/] Tortuosity

%% Verification
    fprintf('The valid range for the Porous material is %g-%g Hz \n',0.01 * sigma, sigma) 
    U = rho_0 * f / sigma;

%% Miki Model
    Z_MIK = rho_0*c_0 * ( 1 + 0.0785*U.^(-0.632) - 1i*0.120*U.^(-0.632) ); 
    km_MIK = omega./c_0 .* ( 1 + 0.122*U.^(-0.618) - 1i * 0.180*U.^(-0.618) );
%% Jac Model    
[Zcf,kf] = ChampouxA1j_coef(omega,phi,sigma,tortu,lambda,lambdap);

%%  Calculation for a rigid Wall
    Z_0=rho_0*c_0; 
  % JAC
    Z_p_jac = -1i.* Zcf .*cot(kf*h);
    alpha_jac = 1 - ( abs( (Zcf-Z_0)./ (Zcf+Z_0) ) ).^2;
    % Miki
    Z_p_MIK = -1i.*Z_MIK .*cot(km_MIK*h);
    alpha_MIK = 1 - ( abs( (Z_p_MIK-Z_0)./( Z_p_MIK +Z_0) ) ).^2;
    
%% Plot(

    if sigma > f(end)
        sigma = 10000;
    end
    figure(1)
    subplot(211)
        plot(f, real([Z_p_jac Z_p_MIK])/(rho_0*c_0))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Re(Z)')
        legend('JAC','MIK')
    
    subplot(212)
        plot(f, imag([Z_p_jac Z_p_MIK])/(rho_0*c_0))
        xlim([0.01 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Imag(Z)')
        legend('JAC','MIK')
    
    figure(2)
        semilogx(f, ([alpha_MIK alpha_jac ]))
        xlim([0.0 * sigma, sigma])
        xlabel('Freq [Hz]')
        ylabel('Absorption ')
        legend('Miki','JAC')
