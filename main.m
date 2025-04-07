clear; clc; close all;
%givens:
Kg = 33.3; %    total gear ratio
Km = 0.0401; %   Motor Constant (V/(rad/sec) or Nm/amp)
J = 0.0005 + (0.2 * (0.2794.^2)) + 0.0015; %    Total Inertia (J = Jhub + J_load + J_extra) (Kgm^2)
Rm = 19.2; %    Armature Resistance (Ohms)
% K1 = 30;
% K3 = 1.5;
K1 = [10 20 5 10 10 10]; %  --gains with changable test values
K3 = [1.5 0 0 1 -1 -0.5]; %CHANGE THIS AND K3 TO GET DIFFERENT PLOTS
s = 0.5;%maybe delete later dont know if this is an unknown

theta_L = 0; % starting angle
theta_D = 0.5; % angle desired


% numerator and denomenator assigning
for i=1:length(K1)
    n1(i) = K1(i)*Kg*Km / (J*Rm); %numerator
    d2(i) = s.^2; %denominator ^2
    d1(i) = s.^1 * ( ((Kg.^2)*(Km.^2) / (J*Rm))  +  (K3(i)*Kg*Km / (J*Rm)) ); % denominator ^1
    d0(i) = K1(i).*Kg*Km / (J*Rm); %denominator ^0

    %V_in
    V_in(i) = K1(i)*(theta_D-theta_L)-s*K3(i)*theta_L;
end

figure;
hold on;

%closed loop system
for i = 1:length(K1)
    num = n1(i);
    den = [d2(i) d1(i) d0(i)];
    sysTF = tf(num, den);
    disp(i)
    [x,t] = step(sysTF);
    plot(t,x);
end
legend("10 0", ...
       "20 0", ...
       "5 0", ...
       "10 1", ...
       "10 -1", ...
       "10 -0.5","Location","best")
xlabel("t")
ylabel("x")

%add natural damping
B = 0;

%todo matlab code:
%   determine gains to meet performance objectives
%   modify gains(if necessary) to meet objectives
%   add in natural damping (B = ?) and compare
