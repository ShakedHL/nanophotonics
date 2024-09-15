%סעיף א 

% Constants
n_air = 1.0;    % Refractive index of air
n1 = 1.5;       % Refractive index of first layer (closest to air)
n2 = 1.6;       % Refractive index of second layer (closest to substrate)
n_substrate = 1.5;  % Refractive index of substrate
wavelength = 1.5e-6;  % Wavelength in meters (1.5 microns)

% Function to calculate reflection coefficient for N pairs
N_values = 1:20;  % Number of layer pairs

reflection_intensity = zeros(1, length(N_values)); % Preallocate array

for i = 1:length(N_values)
    N = N_values(i);
    
    % Calculate individual reflection coefficients
    r_air_to_n1 = (n1 - n_air) / (n1 + n_air);
    r_n1_to_n2 = (n2 - n1) / (n2 + n1);
    r_n2_to_substrate = (n_substrate - n2) / (n_substrate + n2);
    
    % Total reflection amplitude for N layers
    total_r = r_air_to_n1;
    
    for j = 1:N
        total_r = total_r + r_n1_to_n2 + r_n2_to_substrate;
    end
    
    reflection_intensity(i) = total_r^2;  % Intensity reflection coefficient (R = |r|^2)
end

% Plot reflection intensity vs N
figure;
plot(N_values, reflection_intensity, 'o-');
xlabel('Number of Layer Pairs (N)');
ylabel('Reflection Intensity');
title('Reflection Intensity vs Number of Layer Pairs');
grid on;

% Curve fitting to an exponential function
fit_func = @(b, N) b(1) * exp(b(2) * N) + b(3);  % Exponential function for fitting
b0 = [0, 0, 0];  % Initial guess for the fitting parameters
b_fit = nlinfit(N_values, reflection_intensity, fit_func, b0);  % Perform the fitting

% Generate fitted data using the obtained parameters
fitted_reflection = fit_func(b_fit, N_values);

% Plot fitted curve on the same graph
hold on;
plot(N_values, fitted_reflection, 'r-');
legend('Original Data', 'Fitted Curve');
hold off;

% Display the fitting parameters
disp('Fitting parameters (a * exp(b * N) + c):');
disp(b_fit);

%סעיף ב

% Constants
n_air = 1.0;        % Refractive index of air
n1 = 1.5;          % Refractive index of first layer (closest to air)
n2 = 1.6;          % Refractive index of second layer (closest to substrate)
n_substrate = 1.5; % Refractive index of substrate
N = 50;            % Number of layer pairs
wavelength_range = linspace(1e-6, 2e-6, 500); % Wavelength range from 1 to 2 microns

% Angles of incidence (in degrees)
angles = [0, 15, 30, 45, 60];

% Initialize reflection intensity results
reflection_results = zeros(length(angles), length(wavelength_range));

% Calculate reflection intensity for each wavelength and angle
for a = 1:length(angles)
    theta_deg = angles(a);
    theta_rad = deg2rad(theta_deg);
    
    for w = 1:length(wavelength_range)
        wavelength = wavelength_range(w);
        
        % Initialize total reflection coefficients
        r_total_parallel = 0;
        r_total_perpendicular = 0;
        
        for j = 1:N
            % Calculate Fresnel coefficients for the current layer
            [r_parallel, r_perpendicular] = fresnel_coefficients(n1, n2, theta_rad);
            
            % Sum the reflection coefficients over all layers
            r_total_parallel = r_total_parallel + r_parallel;
            r_total_perpendicular = r_total_perpendicular + r_perpendicular;
        end
        
        % Calculate total reflection intensity (average of both polarizations)
        reflection_parallel = abs(r_total_parallel)^2;
        reflection_perpendicular = abs(r_total_perpendicular)^2;
        reflection_results(a, w) = (reflection_parallel + reflection_perpendicular) / 2;
    end
end

% Plot reflection intensity vs wavelength for different angles
figure;
hold on;
for a = 1:length(angles)
    plot(wavelength_range * 1e6, reflection_results(a, :), 'DisplayName', [num2str(angles(a)), '°']);
end
xlabel('Wavelength (\mum)');
ylabel('Reflection Intensity');
title('Reflection Intensity vs Wavelength for Various Angles of Incidence');
legend('show');
grid on;
hold off;

%סעיף ג

% Constants
n_air = 1.0;        % Refractive index of air
n1 = 1.5;          % Refractive index of first layer (closest to air)
n2 = 1.6;          % Refractive index of second layer (closest to substrate)
n_substrate = 1.5; % Refractive index of substrate
N = 50;            % Number of layer pairs
wavelength_resonance = 1.5e-6;  % Desired resonance wavelength (1.5 microns)
wavelength_range = linspace(1e-6, 2e-6, 500); % Wavelength range from 1 to 2 microns

% Angles of incidence (in degrees)
angle_of_incidence = 0;  % Normal incidence

% Calculate thicknesses
d1 = wavelength_resonance / (4 * n1);  % Normal thickness of layer 1 (n = 1.5)
d2 = wavelength_resonance / (4 * n2);  % Normal thickness of layer 2 (n = 1.6)

% Introduce perturbation: increase thickness of the (N/2)th layer
perturbation_index = floor(N/2);  % Perturb the middle layer
perturbation_factor = 1.1;        % 10% increase in thickness
perturbed_thickness = d1 * perturbation_factor;  % Increased thickness for resonance

% Initialize reflection intensity results
reflection_intensity = zeros(1, length(wavelength_range));

% Calculate reflection intensity for each wavelength
for w = 1:length(wavelength_range)
    wavelength = wavelength_range(w);
    
    % Initialize total reflection coefficients
    r_total_parallel = 0;
    r_total_perpendicular = 0;
    
    for j = 1:N
        % Determine the thickness for the current layer
        if j == perturbation_index
            layer_thickness = perturbed_thickness;  % Use the perturbed thickness
        else
            if mod(j, 2) == 1
                layer_thickness = wavelength / (4 * n1);  % Layer 1 (n = 1.5)
            else
                layer_thickness = wavelength / (4 * n2);  % Layer 2 (n = 1.6)
            end
        end
        
        % Calculate Fresnel coefficients for the current layer
        [r_parallel, r_perpendicular] = fresnel_coefficients(n1, n2, deg2rad(angle_of_incidence));
        
        % Sum the reflection coefficients over all layers
        r_total_parallel = r_total_parallel + r_parallel;
        r_total_perpendicular = r_total_perpendicular + r_perpendicular;
    end
    
    % Calculate total reflection intensity (average of both polarizations)
    reflection_parallel = abs(r_total_parallel)^2;
    reflection_perpendicular = abs(r_total_perpendicular)^2;
    reflection_intensity(w) = (reflection_parallel + reflection_perpendicular) / 2;
end

% Plot reflection intensity vs wavelength
figure;
plot(wavelength_range * 1e6, reflection_intensity);
xlabel('Wavelength (\mum)');
ylabel('Reflection Intensity');
title('Reflection Intensity vs Wavelength with Resonance Perturbation');
grid on;

% Function to calculate Fresnel coefficients for parallel and perpendicular polarization
function [r_parallel, r_perpendicular] = fresnel_coefficients(n1, n2, theta1)
    theta2 = asin((n1 / n2) * sin(theta1));
    r_parallel = (n1 * cos(theta1) - n2 * cos(theta2)) / (n1 * cos(theta1) + n2 * cos(theta2));
    r_perpendicular = (n2 * cos(theta1) - n1 * cos(theta2)) / (n2 * cos(theta1) + n1 * cos(theta2));
end
