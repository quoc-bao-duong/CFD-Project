clear;
clc;

% Case: Delta_x=Delta_y

% INPUT
    delta_x=1/80;
    delta_y=delta_x;
    Gamma=1; % Diffusion coefficient
    rho=1;     % Density
% DISCRETIZATION
    [N]=Discretization(delta_x,delta_y);
% BOUNDARY CONDITIONS
    [phi_L phi_T]=BCs(delta_x,N);
% VELOCITY FIELD
    [u v]=Velocity(delta_x,N);
% ---------------------------SOLVER--------------------------------%
% UPWIND 1ST ORDER
    [A,B]=Upwind_1st(N,rho,delta_x,delta_y,Gamma,u,v,phi_T,phi_L);
% CENTRAL
    [A1,B1]=Central(N,rho,delta_x,delta_y,Gamma,u,v,phi_T,phi_L);
% Solver PHI 
    PHI_upwind=Gauss(A,B); 
    PHI_central=Gauss(A1,B1); 
% Result size N x N
    PHI_FIELD_upwind=Result(PHI_upwind,N,phi_T,phi_L);
    PHI_FIELD_central=Result(PHI_central,N,phi_T,phi_L);

% ---------------------------DATA FOR TECPLOT--------------------------------%
    Data_plot(PHI_FIELD_upwind,delta_x,1);  % UPWIND 1ST ORDER
    Data_plot(PHI_FIELD_central,delta_x,2); % CENTRAL

