clc, clearvars
%{
% Parameters
n1 = 1.5;
n2 = 1.6;
lambda = 1.5;
N = 10;

% Define the thickness of each layer based on quarter-wave thickness
thickness1 = lambda / (4 * n1);
thickness2 = lambda / (4 * n2);

% Define the transfer matrix for one period of the stack
M = @(n, thickness) [cos(2*pi*n*thickness/lambda), 1i*sin(2*pi*n*thickness/lambda)/(n);
                     1i*sin(2*pi*n*thickness/lambda)*n, cos(2*pi*n*thickness/lambda)];

% Total transfer matrix
M_total = eye(2);
for i = 1:N
    M_total = M_total * M(n1, thickness1) * M(n2, thickness2);
end

% Calculate reflection coefficient
r = (M_total(1,2) / M_total(1,1));

% Reflection intensity
R = abs(r)^2;

% Plot reflection as function of number of layers N
Ns = 1:50;
R_values = zeros(size(Ns));
for i = 1:length(Ns)
    M_total = eye(2);
    for j = 1:Ns(i)
        M_total = M_total * M(n1, thickness1) * M(n2, thickness2);
    end
    r = (M_total(1,2) / M_total(1,1));
    R_values(i) = abs(r)^2;
end

% Plotting the results
figure;
plot(Ns, R_values);
xlabel('Number of layer pairs (N)');
ylabel('Reflection Intensity');
title('Reflection Intensity vs. Number of Layers (N)');
grid on;

% Fit to a function
p = polyfit(Ns, R_values, 2);
fitted_values = polyval(p, Ns);
hold on;
plot(Ns, fitted_values, 'r--');
legend('Calculated', 'Fitted Curve');

% Fit to a quadratic function
p = polyfit(Ns, R_values, 2);
a = p(1);
b = p(2);
c = p(3);

% Display the fitted coefficients
disp(['Fitted Function: R(N) = ', num2str(a), '*N^2 + ', num2str(b), '*N + ', num2str(c)]);
%}


 %Q1

 clc; clearvars;
clear;
close all;

% Set up parameters
lambda_design = 1.5e-6; % Wavelength in meters
index_high = 1.5; % Refractive index of high-index layer
index_low = 1.6; % Refractive index of low-index layer (adjacent to air)
index_air = 1.0; % Refractive index of air
num_layers = 1:400;
angle_of_incidence = 0;

% Part A: Reflectance vs. Number of Layer Pairs
phase_matrix = [1i, 0; 0, -1i];
interface_air = [1, 1; index_air, -index_air];
inv_interface_air = inv(interface_air);
interface_high = [1, 1; index_high * cos(angle_of_incidence), -index_high * cos(angle_of_incidence)];
inv_interface_high = inv(interface_high);
interface_low = [1, 1; index_low * cos(angle_of_incidence), -index_low * cos(angle_of_incidence)];
inv_interface_low = inv(interface_low);

for i = num_layers
    M = inv_interface_air * (interface_high * phase_matrix * inv_interface_high * interface_low * phase_matrix * inv_interface_low)^i * interface_high;
    reflectance_values(i) = abs(M(2, 1) / M(1, 1))^2;
end

figure(1);
hold on;
x_values = tanh(num_layers * log(index_high / index_low)).^2;
plot(num_layers, reflectance_values);
plot(num_layers, x_values);
title('Reflectance as a function of number of layer pairs');
xlabel('Number of Layer Pairs (N)');
ylabel('Reflectance');
legend('Reflectance', 'tanh(Nx)^2');
grid on;

% Part B: Reflectance for Different Angles of Incidence
num_layers = 50;
angles = [0:15:60];
theta_rad = deg2rad(angles);
wavelength_range = linspace(1e-6, 2e-6, 500);

for j = 1:length(theta_rad)
    reflectance_values = zeros(size(wavelength_range)); % Initialize new matrix for each iteration
    
    theta1 = asin((index_air / index_high) * sin(theta_rad(j)));
    interface_air = [1, 1; index_air * cos(theta_rad(j)), -index_air * cos(theta_rad(j))];
    inv_interface_air = inv(interface_air);
    interface_high = [1, 1; index_high * cos(theta1), -index_high * cos(theta1)];
    inv_interface_high = inv(interface_high);
    theta2 = asin((index_high / index_low) * sin(theta1));
    interface_low = [1, 1; index_low * cos(theta2), -index_low * cos(theta2)];
    inv_interface_low = inv(interface_low);
    theta3 = asin((index_low / index_high) * sin(theta2));
    interface_high_2 = [1, 1; index_high * cos(theta3), -index_high * cos(theta3)];
    inv_interface_high_2 = inv(interface_high_2);
    
    for k = 1:length(wavelength_range)
        exp1 = (lambda_design / wavelength_range(k)) * (pi / 2) * cos(theta1);
        exp2 = (lambda_design / wavelength_range(k)) * (pi / 2) * cos(theta2);
        exp3 = (lambda_design / wavelength_range(k)) * (pi / 2) * cos(theta3);
        
        phase_high = [exp(1i * exp1), 0; 0, exp(-1i * exp1)];
        phase_low = [exp(1i * exp2), 0; 0, exp(-1i * exp2)];
        phase_high_2 = [exp(1i * exp3), 0; 0, exp(-1i * exp3)];
        
        M = inv_interface_air * interface_high * phase_high * inv_interface_high * (interface_low * phase_low * inv_interface_low * interface_high_2 * phase_high_2 * inv_interface_high_2)^(num_layers - 1) * interface_low;
        
        reflectance_values(k) = abs(M(2, 1) / M(1, 1))^2;
    end
    
    figure(2);
    subplot(2, 3, j);
    plot(wavelength_range * 1e6, reflectance_values, 'DisplayName', sprintf('%d°', angles(j)));
    title(sprintf('Reflectance vs. Wavelength for angle = %d°', angles(j)));
    xlabel('Wavelength (\mum)');
    ylabel('Reflectance');
    grid on;
end

% Part C: Reflectance with varying middle layer thickness
clear;
lambda_design = 1.5e-6; % Design wavelength in meters
index_high = 1.5; % High-index layer
index_low = 1.6; % Low-index layer (next to air)
index_air = 1.0; % Air refractive index
num_layers = 50;
angle_of_incidence = 0;
thickness_variation = (0:0.01:2) * 10e-6;

% Initialize matrices
interface_air = [1, 1; index_air, -index_air];
inv_interface_air = inv(interface_air);
interface_high = [1, 1; index_high * cos(angle_of_incidence), -index_high * cos(angle_of_incidence)];
inv_interface_high = inv(interface_high);
interface_low = [1, 1; index_low * cos(angle_of_incidence), -index_low * cos(angle_of_incidence)];
inv_interface_low = inv(interface_low);

reflectance_values = zeros(1, length(thickness_variation)); % Preallocate for reflectance
intensity = zeros(1, length(thickness_variation)); % Preallocate for intensity
transmission = zeros(1, length(thickness_variation)); % Preallocate for transmission

j = 0; % Initialize the index
for k = thickness_variation
    j = j + 1; % Increment the index
    exp1 = (pi / 2);
    exp2 = (2 * pi * index_high / lambda_design) * k;
    
    % Calculate phase matrices
    phase_matrix = [exp(1i * exp1), 0; 0, exp(-1i * exp1)];
    phase_matrix_thick = [exp(1i * exp2), 0; 0, exp(-1i * exp2)];
    
    % Calculate combined matrix
    combined_matrix = interface_high * phase_matrix_thick * inv_interface_high * interface_low * phase_matrix * inv_interface_low;
    M = inv_interface_air * (interface_high * phase_matrix * inv_interface_high * interface_low * phase_matrix * inv_interface_low)^((num_layers / 2) - 1) * combined_matrix * (interface_high * phase_matrix * inv_interface_high * interface_low * phase_matrix * inv_interface_low)^(num_layers / 2) * interface_high;
    
    % Calculate reflectance, transmission, and intensity
    reflectance_values(j) = abs(M(2, 1) / M(1, 1))^2;
    transmission(j) = abs(1 / M(1, 1))^2;
    intensity(j) = (1 / (1 - reflectance_values(j))^2) * (transmission(j)^2 / (1 + (4 * reflectance_values(j) / (1 - reflectance_values(j))^2) * (sin(2 * pi * k * 1e6))^2));
end

% Plot reflectance
figure(3);
hold on;
plot(thickness_variation, reflectance_values);
title('Reflection as a function of middle layer thickness');
xlabel('Thickness of the layer (m)');
ylabel('Reflectance');
legend('Reflectance');
grid on;

% Plot intensity
figure(4);
hold on;
plot(thickness_variation, intensity);
title('Intensity as a function of middle layer thickness');
xlabel('Thickness of the layer (m)');
ylabel('Intensity');
legend('Intensity');
grid on;

% Part D: Reflectance for Different Wavelengths with Thickness Variation
clear;
lambda_design = 1.5e-6;
index_high = 1.5;
index_low = 1.6;
index_air = 1.0;
num_layers = 10;
angle_of_incidence = 0;
wavelength_range = linspace(1e-6, 2e-6, 2000);
layer_thickness = 20e-6;

interface_air = [1, 1; index_air, -index_air];
inv_interface_air = inv(interface_air);
interface_high = [1, 1; index_high * cos(angle_of_incidence), -index_high * cos(angle_of_incidence)];
inv_interface_high = inv(interface_high);
interface_low = [1, 1; index_low * cos(angle_of_incidence), -index_low * cos(angle_of_incidence)];
inv_interface_low = inv(interface_low);

for k = 1:length(wavelength_range)
    exp1 = (lambda_design / wavelength_range(k)) * (pi / 2);
    exp2 = layer_thickness * (2 * pi * index_high / wavelength_range(k));
    
    phase_matrix = [exp(1i * exp1), 0; 0, exp(-1i * exp1)];
    phase_matrix_thick = [exp(1i * exp2), 0; 0, exp(-1i * exp2)];
    
    combined_matrix = interface_high * phase_matrix_thick * inv_interface_high * interface_low * phase_matrix * inv_interface_low;
    M = inv_interface_air * (interface_high * phase_matrix * inv_interface_high * interface_low * phase_matrix * inv_interface_low)^((num_layers / 2) - 1) * combined_matrix * (interface_high * phase_matrix * inv_interface_high * interface_low * phase_matrix * inv_interface_low)^(num_layers / 2) * interface_high;
    
    reflectance_values(k) = abs(M(2, 1) / M(1, 1))^2;
end

figure(5);
plot(wavelength_range, reflectance_values);
title(sprintf('Reflection as a function of middle layer expansion | %d layers', num_layers));
xlabel('Wavelength (m)');
ylabel('Reflectance');
legend('Reflectance');
grid on;


