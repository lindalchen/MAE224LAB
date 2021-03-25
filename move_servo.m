%% Lab 03 
% MAE 224
% Last edited: lindachen@princeton.edu
clear all
close all
clc
%
name = 'mae224-2021-B';
%Enter the unique access token for your photon%
atoken = 'abb13444463d39ffdde80df8690effb884a756d2';
g = Photon(name, atoken);

g.attachServo('D0');

%%
servo_position_ticks = 60; % from 60 to 115 ticks on servo

g.move(servo_position_ticks); 

% lowest positon is at 115 ticks
%% Calibration
tics = [60 70 80 85 90 110]; % ticks
y = [29.5 39.5 49 54.5 60 81]; % mm (distance above the channel surface)
n = 1; % number of degrees 
% polynomial fit
c = polyfit(tics, y, n);

slope = c(1); % mm/ticks

% equation
% h = (115-x)*slope+1
tic_positions_to_test = [60 70 80 85 90 110]; % mm
%% Voltages
% Create a 2D array called "voltages" to store the recorded voltages in. The 
% rows should correspond to the different tick positions to test; the 
% columns should correspond to each of the ~20 measurements you make at
% each position.
rows = length(tic_positions_to_test);
num_measure = 10;
voltages = zeros(rows, num_measure);

% only change these few lines if you're re-starting measurements after the
% code crashed part-way through the data acquisition
starting_tick_i = 1;
if false
    load(data_file_name); % load the tick positions and voltages as previously saved
    starting_tick_i = 15; % change to whatever row in "voltages" you want to start with
end

% Loop through each position x
for i = 1:rows % ... something involving the tick positions to test
    
    % move the servo to the correct position, wait a few seconds
    g.move((tic_positions_to_test(i)));
    
    % loop through the multiple measurements made at this tick position
    for j = 1:num_measure % ... something involving the number of measurements at each position
        
        voltages(i,j) = g.analogRead('A0'); % read the voltage, store it in the correct spot in the
        % "voltages" matrix

        % wait half a second
        pause(0.5);
        
        % if you want, you can also have it calculate and print the 
        % velocity corresponding to the voltage that was just read
        fprintf('\n voltage %.2f \n', voltages(i,j));
        fprintf('pressure %.2f \n', calcurve(voltages(i,j)));
        fprintf('velocity %.2f \n',bernoulli(calcurve(voltages(i,j))));
    end
    
    % save the tick_positions_to_test vector and the voltages array to a
    % file
    
end
save('lab_3_data.mat','tic_positions_to_test','voltages');

function p = calcurve(Vs)
    V_plus = 3.1; % volts
    Pmax = 248.84; % Pa
    Pmin = -248.84; % Pa
    p = Pmin + (Pmax-Pmin)/0.8*((Vs/V_plus)-0.1);
end

function v = bernoulli(p)
    rho = 1.23; % kg/m^3
    v = sqrt(2*p/rho); % m/s
end
