clear;
clc;

% Case: Delta_x=Delta_y

% INPUT
    delta_x=1/80; % ~ delta_x = 0.0714
    delta_y=delta_x;
    delta_t=10^-5;
    t_max=50;
    Gamma=0.5; % Diffusion coefficient
    rho=1;     % Density
    value=0.5;   % Phi initial
% DISCRETIZATION
    [N]=Discretization(delta_x,delta_y);
% BOUNDARY CONDITIONS
    [phi_L phi_T]=BCs(delta_x,N);
% VELOCITY FIELD
    [u v]=Velocity(delta_x,N);
% ---------------------------SOLVER--------------------------------%
% PHI INITIAL ( PHI(X,Y) AT T = 0)
    [Phi_i]=Phi_Initial(value,N);
% FORWARD EULER -EXPLICIT
    [Phi_new]=Forward_Euler(N,rho,delta_x,delta_y,Gamma,u,v,phi_T,phi_L,delta_t,t_max,Phi_i)

% ---------------------------DATA FOR TECPLOT--------------------------------%
    Data_plot(Phi_i,delta_x,0);   % PHI INITIAL 
    Data_plot(Phi_new,delta_x,2); % PHI (X,Y) T=...


