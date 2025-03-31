clear; clc; close all;

Kg = 33.3;
Km = 0.0401; % N*m/amp
J = 0.0005 + 0.2*0.2794^2 + 0.0015; % kg*m^2
Rm = 19.2; % ohms
K1 = 10;
K3 = 0;

n1 = K1*Kg*Km/(J*Rm);
d2 = 1;
d1 = ((Kg*Km)^2 + K3*Kg*Km)/(J*Rm);
d0 = n1;


%% Closed Loop System
num = n1;
den = [d2 d1 d0];
sysTF = tf(num,den);

%% Step Response
[x,t] = step(sysTF);