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

%% Calibration with our own data
servo_position_ticks = 60; % from 60 to 115 ticks on servo

g.move(servo_position_ticks); 
% lowest positon is at 115 ticks
%% Recording the Voltages
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
        fprintf('velocity %.2f \n', bernoulli(calcurve(voltages(i,j))));
    end
    
    % save the tick_positions_to_test vector and the voltages array to a
    % file
    
end
save('lab_3_data.mat','tic_positions_to_test','voltages');

%% Processing with provided data 03/25
% loading
% knob 4 close
num = '4';
pos = 'close';
load(append('speed_', num, '_position_', pos, '.mat'));
close_4_tick_pos = tick_positions_to_test;
close_4_voltages = voltages;

% knob 7 close
num = '7';
load(append('speed_', num, '_position_', pos, '.mat'));
close_7_tick_pos = tick_positions_to_test;
close_7_voltages = voltages;

% knob 10 close
num = '10';
load(append('speed_', num, '_position_', pos, '.mat'));
close_10_tick_pos = tick_positions_to_test;
close_10_voltages = voltages;

% knob 4 far
num = '4';
pos = 'far';
load(append('speed_', num, '_position_', pos, '.mat'));
far_4_tick_pos = tick_positions_to_test;
far_4_voltages = voltages;

% knob 7 far
num = '7';
load(append('speed_', num, '_position_', pos, '.mat'));
far_7_tick_pos = tick_positions_to_test;
far_7_voltages = voltages;

% knob 10 far
num = '10';
load(append('speed_', num, '_position_', pos, '.mat'));
far_10_tick_pos = tick_positions_to_test;
far_10_voltages = voltages;

%%
% Done with calibration based on provided data (see GitHub manual)
slope = 1; % mm/ticks
lowest_tick = 118;
close_4_heights_to_test = (lowest_tick-close_4_tick_pos)*slope; % mm
close_7_heights_to_test = (lowest_tick-close_7_tick_pos)*slope; % mm
close_10_heights_to_test = (lowest_tick-close_10_tick_pos)*slope; % mm

far_4_heights_to_test = (lowest_tick-far_4_tick_pos)*slope; % mm
far_7_heights_to_test = (lowest_tick-far_7_tick_pos)*slope; % mm
far_10_heights_to_test = (lowest_tick-far_10_tick_pos)*slope; % mm

num_measure = 20;

% Velocities
close_4_velocity = bernoulli(calcurve(close_4_voltages)); 
close_7_velocity = bernoulli(calcurve(close_7_voltages));
close_10_velocity = bernoulli(calcurve(close_10_voltages)); 
            
far_4_velocity = bernoulli(calcurve(far_4_voltages)); 
far_7_velocity = bernoulli(calcurve(far_7_voltages)); 
far_10_velocity = bernoulli(calcurve(far_10_voltages)); 

% Average and standard deviation across all rows
close_4_mean_velocity =  mean(close_4_velocity, 2);
close_4_std_velocity = std(close_4_velocity,0, 2);

close_7_mean_velocity =  mean(close_7_velocity, 2);
close_7_std_velocity = std(close_7_velocity,0, 2);

close_10_mean_velocity =  mean(close_10_velocity, 2);
close_10_std_velocity = std(close_10_velocity,0, 2);

far_4_mean_velocity =  mean(far_4_velocity, 2);
far_4_std_velocity = std(far_4_velocity,0, 2);

far_7_mean_velocity =  mean(far_7_velocity, 2);
far_7_std_velocity = std(far_7_velocity,0, 2);

far_10_mean_velocity =  mean(far_10_velocity, 2);
far_10_std_velocity = std(far_10_velocity,0, 2);

% Calculating free stream velcocity
close_4_fs_vel = max(close_4_mean_velocity);
close_7_fs_vel = max(close_7_mean_velocity);
close_10_fs_vel = max(close_10_mean_velocity);

far_4_fs_vel = max(far_4_mean_velocity);
far_7_fs_vel = max(far_7_mean_velocity);
far_10_fs_vel = max(far_10_mean_velocity);

% calculating boundary layer thickness fix thi
[~, close_4_bd_index] = min(abs(close_4_mean_velocity- 0.99*close_4_fs_vel));
[~, close_7_bd_index] = min(abs(close_7_mean_velocity- 0.99*close_7_fs_vel));
[~, close_10_bd_index] = min(abs(close_10_mean_velocity- 0.99*close_10_fs_vel));

close_4_bd = close_4_heights_to_test(close_4_bd_index);
close_7_bd = close_7_heights_to_test(close_7_bd_index);
close_10_bd = close_10_heights_to_test(close_10_bd_index);

[~, far_4_bd_index] = min(abs(far_4_mean_velocity-0.99*far_4_fs_vel));
[~, far_7_bd_index] = min(abs(far_7_mean_velocity-0.99*far_7_fs_vel));
[~, far_10_bd_index] = min(abs(far_10_mean_velocity-0.99*far_10_fs_vel));

far_4_bd = far_4_heights_to_test(far_4_bd_index);
far_7_bd = far_7_heights_to_test(far_7_bd_index);
far_10_bd = far_10_heights_to_test(far_10_bd_index);

bd_height_4 = [close_4_bd far_4_bd];
bd_height_7 = [close_7_bd far_7_bd];
bd_height_10 = [close_10_bd far_10_bd];

x_value = [0.494 1.495]; % m

% Plotting
figure(1)
set(gca, 'FontSize', 16) 
hold on
grid on
errorbar(close_4_mean_velocity, close_4_heights_to_test, close_4_std_velocity);
errorbar(far_4_mean_velocity, far_4_heights_to_test, far_4_std_velocity);

errorbar(close_7_mean_velocity, close_7_heights_to_test, close_7_std_velocity);
errorbar(far_7_mean_velocity, far_7_heights_to_test, far_7_std_velocity);

errorbar(close_10_mean_velocity, close_10_heights_to_test, close_10_std_velocity);
errorbar(far_10_mean_velocity, far_10_heights_to_test, far_10_std_velocity);

legend(sprintf('%.3f m/s', close_4_fs_vel), sprintf('%.3f m/s', far_4_fs_vel), ...
    sprintf('%.3f m/s', close_7_fs_vel), sprintf('%.3f m/s', far_7_fs_vel), ...
    sprintf('%.3f m/s', close_10_fs_vel), sprintf('%.3f m/s', far_10_fs_vel));
ylabel('Height (mm)');
xlabel('Velocity (m/s)');
title('Velocity Profiles')

figure(2)
clf
hold on
grid on
set(gca, 'FontSize', 16)
scatter(x_value, bd_height_4, 'LineWidth', 2);
scatter(x_value, bd_height_7, 'LineWidth', 4);
scatter(x_value, bd_height_10, 'LineWidth', 2);
xlabel('Distance Along Tunnel (m)')
ylabel('Height (mm)');
title('Boundary Layer Thickness of Various Velocities Along Wind Tunnel');
avg_4 = mean([close_4_fs_vel, far_4_fs_vel]);
avg_7 = mean([close_7_fs_vel, far_7_fs_vel]);
avg_10 = mean([close_10_fs_vel, far_10_fs_vel]);
legend(sprintf('%.3f m/s', avg_4), ...
    sprintf('%.3f m/s', avg_7), ...
    sprintf('%.3f m/s', avg_10)); 

%% Discussion
% Non-dimnesionalize velocity profile
close_4_nond_vel = close_4_mean_velocity ./ close_4_fs_vel;
close_7_nond_vel = close_7_mean_velocity ./ close_7_fs_vel;
close_10_nond_vel = close_10_mean_velocity ./ close_10_fs_vel;

far_4_nond_vel = close_4_mean_velocity ./ close_4_fs_vel;
far_7_nond_vel = far_7_mean_velocity ./ far_7_fs_vel;
far_10_nond_vel = far_10_mean_velocity ./ far_10_fs_vel;

% calculation n for each profile
v = 0.00001568; % kinematic viscosity of air, grabbed from website

close_4_n = 0.001*close_4_heights_to_test.*(sqrt(close_4_fs_vel./(v.* 0.494)));
close_7_n = 0.001*close_7_heights_to_test.*(sqrt(close_7_fs_vel./(v.* 0.494)));
close_10_n = 0.001*close_10_heights_to_test.*(sqrt(close_10_fs_vel./(v.* 0.494)));

far_4_n = 0.001*far_4_heights_to_test.*(sqrt(far_4_fs_vel./(v.* 1.495)));
far_7_n = 0.001*far_7_heights_to_test.*(sqrt(far_7_fs_vel./(v.* 1.495)));
far_10_n = 0.001*far_10_heights_to_test.*(sqrt(far_10_fs_vel./(v.* 1.495)));
 
% blasius.txt comparison
load('blasius.txt')
blasius_y = blasius(1:2000);
blasius_x = blasius(2001:4000);
 
figure(3)
hold on
grid on
set(gca, 'FontSize', 14)
plot(blasius_x, blasius_y);
plot(close_4_nond_vel, close_4_n);
plot(close_7_nond_vel, close_7_n);
plot(close_10_nond_vel, close_10_n);
plot(far_4_nond_vel, far_4_n);
plot(far_7_nond_vel, far_7_n);
plot(far_10_nond_vel, far_10_n);
legend('Blasius', sprintf('%.3f m/s', close_4_fs_vel), sprintf('%.3f m/s', close_7_fs_vel), ...
    sprintf('%.3f m/s', close_10_fs_vel), sprintf('%.3f m/s', far_4_fs_vel), ...
    sprintf('%.3f m/s', far_7_fs_vel), sprintf('%.3f m/s', far_10_fs_vel));
xlabel('Non-Dimensional Velocity');
ylabel('$\eta$', 'Interpreter', 'Latex');
title('Non-Dimensionalized Velocity Profiles in Comparison with the Blasius Solution');
 
% experimental vs laminar vs turbulent
x_tube = 0:0.01:1.5;
 
for i = 1:length(x_tube)
    Re_4(i) = (close_4_fs_vel*x_tube(i))/v;
    lam_4_bd(i) = (4.91*x_tube(i))./(Re_4(i).^(0.5));
    turb_4_bd(i) = (0.38*x_tube(i))./(Re_4(i).^(0.2));
end
 
bd_4_exp = [0.011,0.015];
figure(4)
hold on
grid on
set(gca, 'FontSize', 14)
plot(x_tube,lam_4_bd);
plot(x_tube,turb_4_bd);
scatter(x_value, bd_4_exp, 'filled');
legend('Laminar', 'Turbulent', 'Experimental');
xlabel('x (m)')
ylabel('$\delta$(x)', 'Interpreter', 'Latex');
title(sprintf('Boundary Layer Profile for %.3f m/s', close_4_fs_vel))

%% Error Analysis
% uncertainty of the average velocity, pressure terms in it
close_4_pres_err = error_sensor(mean(calcurve(close_4_voltages), 2)); 
close_7_pres_err = error_sensor(mean(calcurve(close_7_voltages), 2));
close_10_pres_err = error_sensor(mean(calcurve(close_10_voltages), 2)); 
            
far_4_pres_err = error_sensor(mean(calcurve(far_4_voltages), 2)); 
far_7_pres_err = error_sensor(mean(calcurve(far_7_voltages), 2)); 
far_10_pres_err = error_sensor(mean(calcurve(far_10_voltages), 2)); 

% Plot
figure(5)
hold on
% Averages and Error Bar
% Close 4
errorbar(close_4_mean_velocity, close_4_heights_to_test, close_4_pres_err, ...
    'horizontal');
scatter_c4v = repelem(close_4_heights_to_test,20);
arr_c4v = reshape(close_4_velocity', 600, 1); % see line 120
scatter(arr_c4v, scatter_c4v, 'filled');

% Far 4
errorbar(far_4_mean_velocity, far_4_heights_to_test, far_4_pres_err, ...
    'horizontal');
scatter_f4v = repelem(far_4_heights_to_test,20);
arr_f4v = reshape(far_4_velocity', 600, 1);
scatter(arr_f4v, scatter_f4v, 'filled');

% Close 7
errorbar(close_7_mean_velocity, close_7_heights_to_test, close_7_pres_err, ...
    'horizontal');
scatter_c7v = repelem(close_4_heights_to_test,20);
arr_c7v = reshape(close_7_velocity', 600, 1);
scatter(arr_c7v, scatter_c7v, 'filled');

% Far 7
errorbar(far_7_mean_velocity, far_7_heights_to_test, far_7_pres_err, ...
    'horizontal');
scatter_f7v = repelem(far_7_heights_to_test,20);
arr_f7v = reshape(far_7_velocity', 38*20, 1);
scatter(arr_f7v, scatter_f7v, 'filled');

% Close 10
errorbar(close_10_mean_velocity, close_10_heights_to_test, close_10_pres_err, ...
    'horizontal');
scatter_c10v = repelem(close_4_heights_to_test,20);
arr_c10v = reshape(close_10_velocity', 600, 1);
scatter(arr_c10v, scatter_c10v, 'filled');

% Far 10
errorbar(far_10_mean_velocity, far_10_heights_to_test, far_10_pres_err, ...
    'horizontal');
scatter_f10v = repelem(far_10_heights_to_test,20);
arr_f10v = reshape(far_10_velocity', 600, 1);
scatter(arr_f10v, scatter_f10v, 'filled');

legend(sprintf('%.3f m/s Average with Error', close_4_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', close_4_fs_vel), ...
    sprintf('%.3f m/s Average with Error', far_4_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', far_4_fs_vel), ...
    sprintf('%.3f m/s Average with Error', close_7_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', close_7_fs_vel), ...
    sprintf('%.3f m/s Average with Error', far_7_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', far_7_fs_vel), ...
    sprintf('%.3f m/s Average with Error', close_10_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', close_10_fs_vel), ...
    sprintf('%.3f m/s Average with Error', far_10_fs_vel), ...
    sprintf('%.3f m/s Raw Velocities', far_10_fs_vel))
ylabel('Height (mm)');
xlabel('Velocity (m/s)');
title('Velocity Profiles')

function sigma_v = error_sensor(p_app)
    p_err_value = 248.84; % pascals, sensor range
    p_err_percent = 0.02; % percentage
    p_err = p_err_value*p_err_percent*2; % pascals

    rho = 1.23; % kg/m^3 density of air
    dv_dp = 1./sqrt(2.*rho.*p_app);
    sigma_v = sqrt(p_err.^2.*(dv_dp).^2);
end

function p = calcurve(Vs)
    V_plus = 3.337; % volts
    Pmax = 248.84; % Pa
    Pmin = -248.84; % Pa
    p = Pmin + (Pmax-Pmin)/0.8*((Vs/V_plus)-0.1);
end

function v = bernoulli(p)
    rho = 1.23; % kg/m^3
    v = sqrt(2*p/rho); % m/s
end
