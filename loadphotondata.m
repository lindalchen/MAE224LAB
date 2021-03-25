function [name,atoken,port] = loadphotondata()
% Loads information of Linda's Photon data for MAE221
% Last edited lindachen@princeton.edu

%Enter your photon name
name = 'mae224-2021-B';
%Enter the unique access token for your photon%
atoken = 'abb13444463d39ffdde80df8690effb884a756d2';
%Set the digital pin for the LED%
port = '/dev/tty.usbmodem144301'; % so it reads directly from usb
end

