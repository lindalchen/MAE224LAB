%% Lab 04
% MAE 224
% Last edited: lindachen@princeton.edu

%Name of the photon inside the airfoil
name = 'MRSAIRFOIL'; 
%Access token
atoken = 'e9cf36b5a026caed895292f11d0d2cb7e201fd12'; 

g = Photon(name,atoken);

[pt,pb,aoa] = g.mrsairfoil; % This command collects the pressure data in pt (pressure top), pb (pressure bottom) and angle of attack (aoa) from a gyro sensor embedded in the airfoil. 
%In real experiments, there is actually a for-loop to execute this command multiple times.

[ xn, yn, u, v ] = airfoil_normals( aoa );

%% Calculating Reynold's Number

% 018 NACA blade

rho_air = 1.2; % kg/m^3 air density
mu_air = 1.8e-5; % Pa*s air dynamic viscosity
u1 = 4.5; % m/s
u2 = 5.5; % m/s
u3 = 8.4; % m/s
L = 30.5*10e-2; % m, chord length

Re_1 = rho_air*u1*L/mu_air;
Re_2 = rho_air*u2*L/mu_air;
Re_3 = rho_air*u3*L/mu_air;

% checked with isa and michelle

% assume incompressible flow --no mach
% start -10 deg to 52 deg
% save as .txt files
% x avis cd to alpha (yaxis)
% save by left clicking
% oppoint view f5, choose angle of attack (operational points)

% how do you find critical angle of attck
% eyeball it off the graph?

%% QBlade Simulation



