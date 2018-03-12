%% LVD Project Body Model Initialization (constant_ks_varying_cs.m)

% Created by Brandon Nguyen, Steven Chung, Sarae Ang, Michael Lam
% Modeling the ground as a constant sine wave input, we can see the variation 
% in displacement for the head, lower body, sprung and unsprung masses over 
% the duration of the simulation

% Will probably want to iterate all the parameters over their realistic 
% values in order to 'tune' them, and see how these varied parameters will
% affect the displacement of the head and lower body. Since we're focusing
% on back pain, I think we will want to focus mainly on the displacement of
% the lower body, while keeping head displacements right below the maximum
% before any damage can happen.

clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot adjustments
set(0,'DefaultLineLineWidth',1.5)
set(0,'DefaultLineColor',[1,1,1])
set(0,'DefaultLineMarkerSize',15)
set(0,'DefaultAxesFontSize',15)
set(0,'DefaultFigureColor',[1,1,1])
set(0,'DefaultTextFontSize',26)
set(0,'DefaultTextInterpreter','latex')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontName','Arial')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Quarter Car Parameters

%for big_count = 1:3
    
%{
    switch(big_count)
        case 1          % First run
            b = 1;
        case 2
            b = 2;
        case 3 
            b = 3.7;
    end
%}
    
    %{
    f1 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    f2 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    f3 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    f4 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    %}
    
    %finput = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    %fft_Zh = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    
    bode_h = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    bode_b = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    bode_s = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    bode_u = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
    
    %% Forcing Function

    StopTime = 10;              % Simulation Stop Time [s]

    t = 0:1:50;                 % Time [s]
    unitstep = t >= 0;          % Unit step
    v = t.*unitstep;            % Vehicle speed (ramp up) [m/s]
    b = 1;                      % Bump length [m]
    lambda = 2*b;               % Wavelength [m]
    f = v/lambda;               % Signal frequency [Hz]
    A = 0.1;                    % Signal amplitude [m]

    omega = 2*pi*f;             % Signal frequency [rad/s]
    
    %% Initial Conditions
    
    m_h = 20;                   % Mass of head [kg] (1-50)
    m_b = 45;                   % Mass of body [kg] (1-100)
    m_s = 240;                  % Sprung mass [kg] (1-500)
    m_u = 36;                   % Unsprung mass [kg] (1-50)
    
    cp = 1360;                  % Thorax and Pelvis damping [N*s/m] (500-2000)
    cc = 1649;                  % Seat cushion damping [N*s/m] (500-2000)
    %cs = 1650;                 % Suspension oil damping [N*s/m] (500-3000, Default @ 1650) 

    kp = 45005;                 % Thorax and Pelvis stiffness [N/m] (2000-50000)
    kc = 20000;                 % Seat cushion stiffness [N/m] (2000-50000)
    ks = 10000;                 % Suspension spring stiffness [N/m] (2000-25000)
    kt = 160000;                % Tire stiffness [N/m] (2000-200000)

    % Tablulate Model Characteristics
    cs_vector = zeros(6,1);
    
    wh_vector = zeros(6,1);
    wb_vector = zeros(6,1);
    ws_vector = zeros(6,1);
    wu_vector = zeros(6,1);
    
    wh_d_vector = zeros(6,1);
    wb_d_vector = zeros(6,1);
    ws_d_vector = zeros(6,1);
    wu_d_vector = zeros(6,1);
    
    cs = 0;         
    i = 0;
    
    while (cs < 3000)
        cs = cs + 500;
        
        i = i + 1;
        cs_vector(i) = cs;  % cs vector

        %% Pre-Processing

        % Resonance Frequencies
        wh = sqrt(kp / m_h) / 2*pi;             % Head natural frequency [Hz]
        wb = sqrt((kp + kc) / m_b) / 2*pi;      % Body natural frequency [Hz]
        ws = sqrt((kc + ks) / m_s) / 2*pi;      % Sprung natural frequency [Hz]
        wu = sqrt((ks + kt) / m_u) / 2*pi;      % Unsprung natural frequency [Hz]

        % Damping Coefficients
        zeta_h = cp / (2*wh*m_h);              
        zeta_b = (cp + cc) / (2*wb*m_b);
        zeta_s = (cc + cs) / (2*ws*m_s);
        zeta_u = cs / (2*wu*m_u);

        % Damped Frequencies
        wh_damped = wh*sqrt(1 - zeta_h^2);      % Damped head natural frequency [Hz]
        wb_damped = wb*sqrt(1 - zeta_b^2);      % Damped body natural frequency [Hz]
        ws_damped = ws*sqrt(1 - zeta_s^2);      % Damped sprung natural frequency [Hz]
        wu_damped = wu*sqrt(1 - zeta_u^2);      % Damped unsprung natural frequency [Hz]
        
        % Tabulated Frequencies in Hz
        wh_vector(i) = wh;
        wb_vector(i) = wb;
        ws_vector(i) = ws;
        wu_vector(i) = wu;

        wh_d_vector(i) = wh_damped;
        wb_d_vector(i) = wb_damped;
        ws_d_vector(i) = ws_damped;
        wu_d_vector(i) = wu_damped;
        
        %% Post Processing
        
        % Simulation and plotting
        % fprintf('Running simulation...\n');
        sim('ref_model');                       % Run simulation for the body model
        
        % Crest Factors 
        extremum_h = max(Zh.data); 
        extremum_b = max(Zb.data);
        extremum_s = max(Zs.data);
        extremum_u = max(Zu.data);

        rms_h = rms(Zh.data);
        rms_b = rms(Zb.data);
        rms_s = rms(Zs.data);
        rms_u = rms(Zu.data);
        rms_y = rms(Zy.data);

        crest_h = extremum_h / rms_h;
        crest_b = extremum_b / rms_b;
        crest_s = extremum_s / rms_s;
        crest_u = extremum_u / rms_u;

        % Ride Comfort
        rms_zy_ddot = rms(Zy_ddot.data);
        rms_zh_ddot = rms(Zh_ddot.data);
        rms_zb_ddot = rms(Zb_ddot.data);
        rms_zs_ddot = rms(Zs_ddot.data);
        rms_zu_ddot = rms(Zu_ddot.data);
        
        % Ride Comfort (Transmissibility X/Y)
        rms_zy_ddot(1) = [];                    % Must empty out or will trigger error
        rms_zh_ddot(1) = [];
        rms_zb_ddot(1) = [];
        zhzy_ddot = rms_zh_ddot./rms_zy_ddot;
        zbzy_ddot = rms_zb_ddot./rms_zy_ddot;
        
        rms_h(1) = [];
        rms_b(1) = [];
        rms_s(1) = [];
        rms_u(1) = [];
        rms_y(1) = [];
        %}
        
        % figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8])

        %{
        % Head
        subplot(2, 2, 1);
        plot(head_disp);
        title('Head Displacement, Zh');
        xlabel('Time [s]');
        ylabel('Displacement [m]');

        % Lower Body
        subplot(2, 2, 2);
        plot(body_disp);
        title('Lower Body Displacement, Zb');
        xlabel('Time [s]');
        ylabel('Displacement [m]');

        % Sprung Mass
        subplot(2, 2, 3);
        plot(sprung_disp);
        title('Sprung Mass Displacement, Zs');
        xlabel('Time [s]');
        ylabel('Displacement [m]');

        % Unsprung Mass
        subplot(2, 2, 4);
        plot(unsprung_disp);
        title('Unsprung Mass Displacement, Zu');
        xlabel('Time [s]');
        ylabel('Displacement [m]');

        % Acceleration Plots
        figure
        plot(Zh_ddot);
        title('Head Acceleration, Zh_ddot');
        xlabel('Time [s]');
        ylabel('Acceleration [m/s^2]');

        figure
        plot(Zb_ddot);
        title('Body Acceleration, Zb_ddot');
        xlabel('Time [s]');
        ylabel('Acceleration [m/s^2]');
        %}

        %% Critical Plots
        
        %{
        % RMS Acceleration vs Vehicle Velocity
        %f1 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
        figure(f1)
        plot(v, rms_zh_ddot);
        title(sprintf('Ride Comfort vs Vehicle Velocity: Bump Length %0.1f m', b));
        xlabel('Vehicle Velocity [m/s]');
        ylabel('RMS Head Acceleration [m/s^2]');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on

        %f2 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
        figure(f2)
        plot(v, rms_zb_ddot);
        title(sprintf('Ride Comfort vs Vehicle Velocity: Bump Length %0.1f m', b));
        xlabel('Vehicle Velocity [m/s]');
        ylabel('RMS Body Acceleration [m/s^2]');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on

        % Displacement vs Vehicle Velocity
        %f3 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
        figure(f3)
        plot(v, rms_h);
        title(sprintf('Position vs Vehicle Velocity: Bump Length %0.1f m', b));
        xlabel('Vehicle Velocity [m/s]');
        ylabel('RMS Head Displacement [m]');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on

        %f4 = figure('WindowStyle','docked','units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
        figure(f4)
        plot(v, rms_b);
        title(sprintf('Position vs Vehicle Velocity:  Bump Length %0.1f m', b));
        xlabel('Vehicle Velocity [m/s]');
        ylabel('RMS Body Displacement [m]');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on
        %}
        
        %% FFT (Transmissibility)
        
        %{
        
        Fs = 200;               % Sampling frequency                    
        T = 1/Fs;               % Sampling period       
        L = length(tout)*2-1;   % Length of signal
        t = (0:L-1)*T;          % Time vector
        f = Fs*(0:(L/2))/L;     % Frequency
        
        Zh = Zh.data;
        
        y = abs(fft(Zh));
        
        figure(fft_Zh)
        plot(f, y);
        title(sprintf('Ride Comfort vs Frequency:  Bump Length %0.1f m', b));
        xlabel('Frequency [Hz]');
        ylabel('Amplitude');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on
        
        %}
        
        %% Data Tabulation
        
        %{
        
        % Print Data
        fprintf('\nNatural frequencies [Hz]: \n');
        fprintf('  Head natural frequency: %0.4f\n', wh);
        fprintf('  Body natural frequency: %0.4f\n', wb);
        fprintf('  Sprung natural frequency: %0.4f\n', ws);
        fprintf('  Unsprung natural frequency: %0.4f\n', wu);

        fprintf('\nDamped frequencies [Hz]: \n');
        fprintf('  Head damped frequency: %0.4f\n', wh_damped);
        fprintf('  Body damped frequency: %0.4f\n', wb_damped);
        fprintf('  Sprung damped frequency: %0.4f\n', ws_damped);
        fprintf('  Unsprung damped frequency: %0.4f\n', wu_damped);

        fprintf('\nCrest Factors: \n');
        fprintf('  Head: %0.4f\n', crest_h);
        fprintf('  Body: %0.4f\n', crest_b);
        fprintf('  Sprung: %0.4f\n', crest_s);
        fprintf('  Unsprung: %0.4f\n', crest_u);

        %}
        
        %% Transfer Function Setup

        %{
        % Case: b = 1 [m] (Includes Road Input Signal)
        num_h = [1.12727*10^11*cs*v (5.09759*10^12*cs*v + 1.12727*10^11*ks*v) (4.5244*10^13*cs*v + 5.09759*10^12*ks*v) 4.5244*10^13*ks*v];
        den_h = [7.776*10^6 (1.1021508*10^9 + 248400*cs)...
                (8.8563*10^10 + 3.4985*10^7*cs + 248400*ks + 7.6746*10^7*v.^2)...
                (6.01197*10^12 + 1.84453*10^9*cs + 3.4985*10^7*ks + 1.08778*10^10*v.^2 + 2.45161*10^6*cs*v.^2)...
                (2.49896*10^14 + 5.40027*10^10*cs + 1.84453*10^9*ks + 8.74082*10^11*v.^2 + 3.45288*10^8*cs*v.^2 + 2.45161*10^6*ks*v.^2)...
                (4.94897*10^15 + 1.19781*10^12*cs + 5.40027*10^10*ks + 5.93357*10^13*v.^2 + 1.82048*10^10*cs*v.^2 + 3.45288*10^8*ks*v.^2)...
                (4.39249*10^16 + 1.62261*10^13*cs + 1.19781*10^12*ks + 2.46638*10^15*v.^2 + 5.32985*10^11*cs*v.^2 + 1.82048*10^10*ks*v.^2)...
                (1.44016*10^14*cs + 1.62261*10^13*ks + 4.88443*10^16*v.^2 + 1.18219*10^13*cs*v.^2 + 5.32985*10^11*ks*v.^2)...
                (1.44016*10^14*ks + 4.33521*10^17*v.^2 + 1.60145*10^14*cs*v.^2 + 1.18219*10^13*ks*v.^2)...
                (1.42138*10^15*cs*v.^2 + 1.60145*10^14*ks*v.^2)...
                1.42138*10^15*ks*v.^2];
        %}
        
        % Case: (Does NOT Include Road Input Signal)
        % Head Transfer Function
        num_h = [3.58822*10^11*cs (1.62261*10^13*cs + 3.58822*10^11*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14*ks];
        den_h = [7.776*10^6 (1.10215*10^9 + 248400*cs)...
                (8.8563*10^10 + 3.4985*10^7*cs + 248400*ks)...
                (6.01197*10^12 + 1.84453*10^9*cs + 3.4985*10^7*ks)...
                (2.49896*10^14 + 5.40027*10^10*cs + 1.84453*10^9*ks)...
                (4.94897*10^15 + 1.19781*10^12*cs + 5.40027*10^10*ks)...
                (4.39249*10^16 + 1.62261*10^13*cs + 1.19781*10^12*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14 ks];
            
        % Body Transfer Function
        num_b = [5.2768*10^9*cs (4.22822*10^11*cs + 5.2768*10^9*ks)...
                (1.62261*10^13*cs + 4.22822*10^11*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14*ks];
        den_b = [7.776*10^6 (1.10215*10^9 + 248400*cs)...
                (8.8563*10^10 + 3.4985*10^7*cs + 248400*ks)...
                (6.01197*10^12 + 1.84453*10^9*cs + 3.4985*10^7*ks)...
                (2.49896*10^14 + 5.40027*10^10*cs + 1.84453*10^9*ks)...
                (4.94897*10^15 + 1.19781*10^12*cs + 5.40027*10^10*ks)...
                (4.39249*10^16 + 1.62261*10^13*cs + 1.19781*10^12*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14*ks];   
            
        % Suspension Transfer Function
        num_s = [1.44*10^8*cs (1.94208*10^10*cs + 1.44*10^8*ks)...
                (8.90874*10^11*cs + 1.94208*10^10*ks)...
                (1.62261*10^13*cs + 8.90874*10^11*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14 ks];
        den_s = [7.776*10^6 (1.10215*10^9 + 248400*cs)...
                (8.8563*10^10 + 3.4985*10^7*cs + 248400*ks)...
                (6.01197*10^12 + 1.84453*10^9*cs + 3.4985*10^7*ks)...
                (2.49896*10^14 + 5.40027*10^10*cs + 1.84453*10^9*ks)...
                (4.94897*10^15 + 1.19781*10^12*cs + 5.40027*10^10*ks)...
                (4.39249*10^16 + 1.62261*10^13*cs + 1.19781*10^12*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                (1.44016*10^14*ks)];

        % Tire Transfer Function
        num_u = [3.456*10^10 (4.89845*10^12 + 1.44*10^8*cs)...
                (2.40013*10^14 + 1.94208*10^10*cs + 1.44*10^8*ks)...
                (4.94897*10^15 + 8.90874*10^11*cs + 1.94208*10^10*ks)...
                (4.39249*10^16 + 1.62261*10^13*cs + 8.90874*10^11*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14*ks];
        den_u = [7.776*10^6 (1.10215*10^9 + 248400*cs)...
                (8.8563*10^10 + 3.4985*10^7*cs + 248400*ks)...
                (6.01197*10^12 + 1.84453*10^9*cs + 3.4985*10^7*ks)...
                (2.49896*10^14 + 5.40027*10^10*cs + 1.84453*10^9*ks)...
                (4.94897*10^15 + 1.19781*10^12*cs + 5.40027*10^10*ks)...
                (4.39249*10^16 + 1.62261*10^13*cs + 1.19781*10^12*ks)...
                (1.44016*10^14*cs + 1.62261*10^13*ks)...
                1.44016*10^14*ks];

        Gh = tf(num_h, den_h);
        Gb = tf(num_b, den_b);
        Gs = tf(num_s, den_s);
        Gu = tf(num_u, den_u);

        figure(bode_h)
        bode(Gh);
        title('Bode Plot of Head Element');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on
        
        figure(bode_b)
        bode(Gb);
        title('Bode Plot of Body Element');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on
        
        figure(bode_s)
        bode(Gs);
        title('Bode Plot of Suspension Element');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on
        
        figure(bode_u)
        bode(Gu);
        title('Bode Plot of Tire Element');
        legend('b1 = 500 [N*s/m]', 'b2 = 1000 [N*s/m]', 'b3 = 1500 [N*s/m]', 'b4 = 2000 [N*s/m]', 'b5 = 2500 [N*s/m]', 'b6 = 3000 [N*s/m]', 'Location', 'best');
        grid on
        hold on

    end
    
    %{
    % Input Signal Information
    figure(finput)
    plot(Zy);
    title(sprintf('Input Signal: Bump Length %0.1f m', b));
    xlabel('Time [s]');
    ylabel('Road Amplitude [m]');
    grid on
    hold on
    %}
    
%end

%{
% Ride Comfort/RMS vs. Damping
figure('units','normalized','outerposition',[0.1 0.1 0.8 0.8]);
plot(cs_vector, rms_zh_ddot);
title(sprintf('Ride Comfort vs Vehicle Velocity: Bump Length %0.1f m', b));
xlabel('Vehicle Velocity [m/s]');
ylabel('RMS Head Acceleration [m/s^2]');
grid on
hold on
%}

%% End

fprintf('Run Successful! \n');