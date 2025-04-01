clear; clc; close all;

%givens:
Kg = 33.3; %    total gear ratio
Km = 0.0401; %   Motor Constant (V/(rad/sec) or Nm/amp)
J = 0.0005 + (0.2 * (0.2794.^2)) + 0.0015; %    Total Inertia (J = Jhub + J_load + J_extra) (Kgm^2)
Rm = 19.2; %    Armature Resistance (Ohms)
K1 = 10; %  --gains with changable test values
K3 = 0;
rad = 0.5; % 

% numerator and denomenator assigning
n1 = K1*Kg*Km / (J*Rm); %numerator
d2 = 1; %denominator ^2
d1 = ((Kg.^2)*(Km.^2) / (J*Rm))  +  (K3*Kg*Km / (J*Rm)); % denominator ^1
d0 = K1*Kg*Km / (J*Rm); %denominator ^0


%closed loop system
num = n1;
den = [d2 d1 d0];
sysTF = tf(num, den);


%Step response
[x,t] = step(sysTF*rad);

figure()
plot(t,x)

%add natural damping
B = 0;

%todo matlab code:
%   determine gains to meet performance objectives
%   modify gains(if necessary) to meet objectives
%   add in natural damping (B = ?) and compare
