clear; clc; close all;
%givens:
Kg = 33.3; %    total gear ratio
Km = 0.0401; %   Motor Constant (V/(rad/sec) or Nm/amp)
J = 0.0005 + (0.2 * (0.2794.^2)) + 0.0015; %    Total Inertia (J = Jhub + J_load + J_extra) (Kgm^2)
Rm = 19.2; %    Armature Resistance (Ohms)
data_import = readMatrix("")

K1 = [10 20 5 10 10 10 15]; %  --gains with changable test values
K3 = [0 0 0 1 -1 -0.5 1]; %CHANGE THIS AND K3 TO GET DIFFERENT PLOTS
s = 0.5;%maybe delete later dont know if this is an unknown
labels = cell(1, length(K1)+1); %prealocating size of labels for legends
theta_L = 0; % starting angle
theta_D = 0.5; % angle desired

s = 0.5; % s for coefficients

T_f = 10;
delta_t = 0.01;
t = (0:delta_t:T_f)';

figure; %plots
hold on;
% numerator and denomenator assigning
for i=1:length(K1)
    n1(i) = K1(i)*Kg*Km / (J*Rm); %numerator
    d2(i) = s.^2; %denominator ^2
    d1(i) = s.^1 * ( ((Kg.^2)*(Km.^2) / (J*Rm))  +  (K3(i)*Kg*Km / (J*Rm)) ); % denominator ^1
    d0(i) = K1(i).*Kg*Km / (J*Rm); %denominator ^0

    %V_in
    V_in(i) = K1(i)*(theta_D-theta_L)-s*K3(i)*theta_L;
end

amplitude = 0.5;% rad, amplitude
T = 10; % s perioid
ref_val  = amplitude * square(2*pi*t/T);

%closed loop system
for i = 1:length(K1)
    sysTF = tf(n1(i), [d2(i) d1(i) d0(i)]); %first is numerator, second is denominator
    %[y_step,t] = step(sysTF);
    y = lsim(sysTF, ref_val, t);
    plot(t,y, 'LineWidth', 1.5);
    labels{i} = sprintf('K1=%g, K3=%g', K1(i), K3(i));
end

plot(t, ref_val,  'k', 'LineWidth', 1);
labels{end} = sprintf('Reference = %g rad', theta_D);
%labels = arrayfun(@(p,d) sprintf('%g, %g', p, d), K1, K3, 'UniformOutput', false);
legend(labels, 'Location','best');
xlabel("Time (sec)")
ylabel("x (rad)")

%add natural damping
% B = 0;


