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

% Generate fitted data using the fitted parameters
fitted_reflection = fit_func(b_fit, N_values);

% Plot fitted curve on the same graph
hold on;
plot(N_values, fitted_reflection, 'r-');
legend('Original Data', 'Fitted Curve');
hold off;

% Display the fitting parameters
disp('Fitting parameters (a * exp(b * N) + c):');
disp(b_fit);

fgszirfhirgrsfkuj
