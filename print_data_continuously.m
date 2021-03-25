%% 
% Last edited: lindachen@princeton.edu

% Setting up the photon
name = 'mae224-2021-B';
%Enter the unique access token for your photon%
atoken = 'abb13444463d39ffdde80df8690effb884a756d2';
g = Photon(name, atoken);

%%
clc
tic;
duration = 100;

while toc < duration
    volt= g.analogRead('A0');
    fprintf('\n voltage %3.2f \n', volt);
    fprintf('pressure %3.2f \n', calcurve(volt));
    fprintf('velocity %3.2f \n',bernoulli(calcurve(volt)));
    pause(2)
    
end

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
