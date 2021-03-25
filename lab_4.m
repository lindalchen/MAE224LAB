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
