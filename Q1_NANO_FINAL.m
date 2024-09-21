clc, clearvars

clear;
close all;

% Define constants
lambda_design = 1.5e-6;    % Wavelength {m}
n_high = 1.5;              
n_low = 1.6;               
n_air = 1.0;               % Refractive index of air
num_pairs = 1:400;         
incident_angle = 0;        

%% Section A: Reflectance vs Number of Layer Pairs
phase_matrix = [1i, 0; 0, -1i];  
D_air = [1, 1; n_air, -n_air];   
D_air_inv = inv(D_air);          
D_high = [1, 1; n_high * cos(incident_angle), -n_high * cos(incident_angle)];
D_high_inv = inv(D_high);
D_low = [1, 1; n_low * cos(incident_angle), -n_low * cos(incident_angle)];
D_low_inv = inv(D_low);

% Loop through pairs to calculate reflectance
for n = num_pairs
    M = D_air_inv * (D_high * phase_matrix * D_high_inv * D_low * phase_matrix * D_low_inv)^n * D_high;
    reflection(n) = abs(M(2, 1) / M(1, 1))^2;
end

% Plotting reflectance as a function of the number of pairs
figure(1);
hold on;
x_values = tanh(num_pairs * log(n_high / n_low)).^2;
plot(num_pairs, reflection);
plot(num_pairs, x_values);
title('Reflectance vs Number of Layer Pairs');
xlabel('Number of Layer Pairs');
ylabel('Reflectance');
legend('Calculated Reflectance', 'Tanh Approximation');
grid on;

%% Section B: Reflectance vs Wavelength for Different Angles
num_pairs = 50;                        
angles_deg = [0:15:60];                 
angles_rad = deg2rad(angles_deg);       
lambda_range = linspace(1e-6, 2e-6, 500); 
for j = 1:length(angles_rad)
    reflectance_wavelength = zeros(size(lambda_range)); 

    theta1 = asin((n_air / n_high) * sin(angles_rad(j)));
    D_air = [1, 1; n_air * cos(angles_rad(j)), -n_air * cos(angles_rad(j))];
    D_air_inv = inv(D_air);
    D_high = [1, 1; n_high * cos(theta1), -n_high * cos(theta1)];
    D_high_inv = inv(D_high);
    theta2 = asin((n_high / n_low) * sin(theta1));
    D_low = [1, 1; n_low * cos(theta2), -n_low * cos(theta2)];
    D_low_inv = inv(D_low);
    theta3 = asin((n_low / n_high) * sin(theta2));
    D_high_2 = [1, 1; n_high * cos(theta3), -n_high * cos(theta3)];
    D_high_2_inv = inv(D_high_2);

    for k = 1:length(lambda_range)
        exp1 = (lambda_design / lambda_range(k)) * (pi / 2) * cos(theta1);
        exp2 = (lambda_design / lambda_range(k)) * (pi / 2) * cos(theta2);
        exp3 = (lambda_design / lambda_range(k)) * (pi / 2) * cos(theta3);

        P1 = [exp(1i * exp1), 0; 0, exp(-1i * exp1)];
        P2 = [exp(1i * exp2), 0; 0, exp(-1i * exp2)];
        P3 = [exp(1i * exp3), 0; 0, exp(-1i * exp3)];

        M = D_air_inv * D_high * P1 * D_high_inv * (D_low * P2 * D_low_inv * D_high_2 * P3 * D_high_2_inv)^(num_pairs-1) * D_low * P2 * D_low_inv * D_high_2;
        reflectance_wavelength(k) = abs(M(2, 1) / M(1, 1))^2;
    end

    % Plotting results for each angle
    figure(2);
    subplot(2, 3, j);
    plot(lambda_range * 1e6, reflectance_wavelength, 'DisplayName', sprintf('%d°', angles_deg(j)));
    title(sprintf('Reflectance vs Wavelength for Angle = %d°', angles_deg(j)));
    xlabel('Wavelength (\mum)');
    ylabel('Reflectance');
    grid on;
end

%% Section C: Reflection with Layer Expansion
clear;
lambda_design = 1.5e-6; 
n_high = 1.5;           
n_low = 1.6;            
n_air = 1.0;            
num_pairs = 50;         
theta_incident = 0;     
layer_thickness = (0:0.01:2) * 10e-6;

D_air = [1, 1; n_air, -n_air];
D_air_inv = inv(D_air);
D_high = [1, 1; n_high * cos(theta_incident), -n_high * cos(theta_incident)];
D_high_inv = inv(D_high);
D_low = [1, 1; n_low * cos(theta_incident), -n_low * cos(theta_incident)];
D_low_inv = inv(D_low);
index = 0;

for thickness = layer_thickness
    index = index + 1;
    exp_high = (pi / 2);
    exp_low = (2 * pi * n_high / lambda_design) * thickness;

    P_high = [exp(1i * exp_high), 0; 0, exp(-1i * exp_high)];
    P_low = [exp(1i * exp_low), 0; 0, exp(-1i * exp_low)];

    M_change = D_high * P_low * D_high_inv * D_low * P_high * D_low_inv;
    M_final = D_air_inv * (D_high * P_high * D_high_inv * D_low * P_high * D_low_inv)^((num_pairs / 2) - 1) * M_change * (D_high * P_high * D_high_inv * D_low * P_high * D_low_inv)^(num_pairs / 2) * D_high;

    reflectance(index) = abs(M_final(2, 1) / M_final(1, 1))^2;

    T(index) = abs(1 / M_final(1, 1))^2;
    a = 1 / (1 - reflectance(index))^2;
    b = 1 + (4 * reflectance(index) / (1 - reflectance(index))^2) * (sin(2 * pi * thickness * 10^6))^2;
    I(index) = a * (T(index)^2 / b);
end

figure(3);
hold on;
plot(layer_thickness, reflectance);
title('Reflectance as a function of middle layer thickness');
xlabel('Layer Thickness');
ylabel('Reflectance');
grid on;

figure(4);
hold on;
plot(layer_thickness, I);
title('I as a function of middle layer thickness');
xlabel('Layer Thickness');
ylabel('I');
grid on;

%% Section D: Reflection vs Wavelength
clear;
lambda_design = 1.5e-6; 
n_high = 1.5;           
n_low = 1.6;            
n_air = 1.0;          
num_pairs = 10;         
theta = 0;              
wavelength_range = linspace(1e-6, 2e-6, 2000);
d_layer = 20 * 10^-6;

D_air = [1, 1; n_air, -n_air];
D_air_inv = inv(D_air);
D_high = [1, 1; n_high * cos(theta), -n_high * cos(theta)];
D_high_inv = inv(D_high);
D_low = [1, 1; n_low * cos(theta), -n_low * cos(theta)];
D_low_inv = inv(D_low);
counter = 0;

for k = 1:length(wavelength_range)
    counter = counter + 1;
    exp1 = (lambda_design / wavelength_range(k)) * (pi / 2);
    exp2 = d_layer * (2 * pi * n_high / wavelength_range(k));

    P1 = [exp(1i * exp1), 0; 0, exp(-1i * exp1)];
    P_change = [exp(1i * exp2), 0; 0, exp(-1i * exp2)];

    M_mod = D_high * P_change * D_high_inv * D_low * P1 * D_low_inv;
    M_reflection = D_air_inv * (D_high * P1 * D_high_inv * D_low * P1 * D_low_inv)^((num_pairs / 2) - 1) * M_mod * (D_high * P1 * D_high_inv * D_low * P1 * D_low_inv)^(num_pairs / 2) * D_high;
    reflectance(counter) = abs(M_reflection(2, 1) / M_reflection(1, 1))^2;
end

figure(5);
plot(wavelength_range, reflectance);
title(sprintf('Reflectance vs Wavelength for %d layers', num_pairs));
xlabel('Wavelength (\mum)');
ylabel('Reflectance');
grid on;
