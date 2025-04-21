    clear; clc; close all;
    %givens:
    Kg = 33.3; %    total gear ratio
    Km = 0.0401; %   Motor Constant (V/(rad/sec) or Nm/amp)
    J = 0.0005 + (0.2 * (0.2794.^2)) + 0.0015; %    Total Inertia (J = Jhub + J_load + J_extra) (Kgm^2)
    Rm = 19.2; %    Armature Resistance (Ohms)
    data_import = readmatrix('data_2.txt');
    time_import = data_import(:,1);
    
    time_import = time_import - time_import(1); %normalizes the time
    idx = (time_import >= 5*1000) & (time_import <= 15*1000);
    
    %PLOT: Plot of hardware data
    figure; hold on;

    title('Hardware x position & Reference vs Time (5–15 sec)')
    plot((time_import(idx)/1000)-5, data_import(idx,2))
    plot((time_import(idx)/1000)-5, data_import(idx,6))
    xlabel('time (s)')
    ylabel('x (rad)')
    legend('actual', 'reference data','Location','best')
    hold off;
    
    
    K1 = [15 10 20 5 10 10 10]; %  --gains with changable test values
    K3 = [1 0 0 0 1 -1 -0.5]; %CHANGE THIS AND K3 TO GET DIFFERENT PLOTS
    s = 0.5;%maybe delete later dont know if this is an unknown
    labels = cell(1, length(K1)+1); %prealocating size of labels for legends
    theta_L = 0; % starting angle
    theta_D = 0.5; % angle desired
    
    s = 0.5; % s for coefficients
    
    T_f = 10;
    delta_t = 0.01;
    t = (0:delta_t:T_f)';
    %PLOT: plot of all functions on same figure;
    figure; 
    hold on;
    % numerator and denomenator assigning
    for i=1:length(K1)
        n1(i) = K1(i)*Kg*Km / (J*Rm); %numerator
        d2(i) = s.^2; %denominator ^2
        d1(i) = s.^1 * ( ((Kg.^2)*(Km.^2) / (J*Rm))  +  (K3(i)*Kg*Km / (J*Rm)) ); % denominator ^1
        d0(i) = K1(i).*Kg*Km / (J*Rm); %denominator ^0
    
        %V_in
        V_in(i) = K1(i)*(theta_D-theta_L)-s*K3(i)*theta_L;%calculates the voltage for each value, makes sure that voltage is not over limit
    end
    
    amplitude = 0.5;% rad, amplitude
    T = 10; % s period for wave
    ref_val  = amplitude*sign(sin(2*pi/T*t));%computes the reference wave using a sine wave where if the val is GEQ 0, = 1, LEQ 0, -1
    
    
    for i = 1:length(K1)
        sysTF = tf(n1(i), [d2(i) d1(i) d0(i)]); %first is numerator, second is denominator
        y = lsim(sysTF, ref_val, t);%key part, makes sysTF attempt to follow ref_val (similar to actual hardware)
        plot(t,y, 'LineWidth', 1.5);
        labels{i} = sprintf('K_1=%g, K_3=%g', K1(i), K3(i));%creates labels for legend
    end
    
    plot(t, ref_val,  'k', 'LineWidth', 1);%plot for reference value
    labels{end} = sprintf('Ref = %g rad', theta_D);%add reference value to legend
    title('Closed-Loop response for K1 and K3 values vs time')
    legend(labels, 'Location','best');
    xlabel("Time (sec)")
    ylabel("x (rad)")
    
    %PLOT: plot of all values on their own figure using subplot
    figure;
    hold on
    sgtitle('Closed-Loop response for K1 and K3 values vs time')
    for i = 1:(length(K1))
        sysTF = tf(n1(i), [d2(i) d1(i) d0(i)]); %first is numerator, second is denominator
        y = lsim(sysTF, ref_val, t);%key part, makes sysTF attempt to follow ref_val (similar to actual hardware)
        subplot(2,4,i)
        title(sprintf('K_1=%g, K_3=%g', K1(i), K3(i)))
        hold on;
        plot(t,y, 'LineWidth', 1.5);
        plot(t, ref_val,  'k', 'LineWidth', 1);%plot for reference value
        %leg = legend(labels{i}, labels{end}, 'Location','best');
        %leg.FontSize = 5; % in case you want a legend with it, not necessary b/c titles cover all info needed   
        xlabel("Time (sec)")
        ylabel("x (rad)")
        hold off;
    end
    hold off;
    
    %PLOT: plot of comparison between computed and Hardware data for K1 = 15. K3 = 1
    figure;
    hold on;
    title('Hardware vs Simulated Output (K_1=15, K_3=1)')

    plot((time_import(idx)/1000)-5, data_import(idx,2))
    plot((time_import(idx)/1000)-5, data_import(idx,6),'LineWidth',1.25)
    
    sysTF = tf(n1(1), [d2(1) d1(1) d0(1)]); %first is numerator, second is denominator (this one resets the previous values working only if the first index is 15 1
    y = lsim(sysTF, ref_val, t);
    plot(t, y, 'LineWidth', 1.5);
    legend('actual', 'reference data', labels{1},'Location','best')
    hold off;