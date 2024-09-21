%   q1  

clc, clearvars

% Parameters
n1 = 1.5;  % Refractive index of layer 1
n2 = 1.6;  % Refractive index of layer 2
lambda = 1.5;  % Wavelength in microns
N = 10;  % Number of layer pairs (can vary this value)

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
p = polyfit(Ns, R_values, 2);  % p will contain [a, b, c]
a = p(1);
b = p(2);
c = p(3);

% Display the fitted coefficients
disp(['Fitted Function: R(N) = ', num2str(a), '*N^2 + ', num2str(b), '*N + ', num2str(c)]);

%Q2 מורחבת 

%{
% Parameters
lambda_min = 500e-9;  % Minimum wavelength in meters (500 nm)
lambda_max = 1500e-9; % Maximum wavelength in meters (1500 nm)
lambda = linspace(lambda_min, lambda_max, 100); % Create a range of wavelengths
% 11 valuse between 500 nm and 1500 nm with steps of 50 nm are :([ 500.,  600.,  700.,  800.,  900., 1000., 1100., 1200., 1300.,
      % 1400., 1500.]) 

      % Define wavelength values in nanometers
wavelengths_nm = [500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500];

% Real part of the dielectric constant (epsilon1) for silver
epsilon_real = [-16.0, -15.5, -14.0, -12.5, -11.5, -10.0, -9.0, -8.5, -8.0, -7.8, -7.5];

% Imaginary part of the dielectric constant (epsilon2) for silver
epsilon_imag = [1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0];

% Combine real and imaginary parts into complex dielectric constant
epsilon_complex = epsilon_real + 1i * epsilon_imag;

% Display the results
disp('Wavelength (nm)    Dielectric Constant (Complex)');
for i = 1:length(wavelengths_nm)
    fprintf('%5.0f              %7.2f + %.2fj\n', wavelengths_nm(i), real(epsilon_complex(i)), imag(epsilon_complex(i)));
end

% Plot the real and imaginary parts of the dielectric constant as a function of wavelength
figure;
subplot(2, 1, 1);
plot(wavelengths_nm, epsilon_real, 'b', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Real Part (\epsilon_1)');
title('Real Part of Dielectric Constant for Silver');
grid on;

subplot(2, 1, 2);
plot(wavelengths_nm, epsilon_imag, 'r', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Imaginary Part (\epsilon_2)');
title('Imaginary Part of Dielectric Constant for Silver');
grid on;

% Constants
c = 3e8;  % Speed of light in vacuum (m/s)
e0 = 8.854e-12; % Permittivity of free space
k_air = 1; % Relative permittivity of air
k_metal = -10 + 1i*0.1; % Relative permittivity of metal (example)

% Calculate the plasmon wavelength
k_plasmon = sqrt((k_air * k_metal) / (k_air + k_metal));  % Dispersion relation

% Propagation distance
propagation_distance = (lambda ./ real(k_plasmon));

% Normalized propagation distance
normalized_distance = propagation_distance ./ lambda;

% Plot results
figure;
plot(lambda*1e9, normalized_distance, 'b', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Normalized Propagation Distance');
title('Normalized Propagation of Surface Plasmons at Metal-Air Interface');
grid on;
hold on;

% Highlight the specific points for lambda = 500 nm and 1500 nm
scatter([500, 1500], [normalized_distance(1), normalized_distance(end)], 'r', 'filled');

% Correct the legend
legend('Normalized Propagation Distance', 'Selected Points: \lambda = 500 nm and 1500 nm');
hold off;

%%ניסיון 2

% Define wavelength values in nanometers (exactly 11 values)
wavelengths_nm = linspace(500, 1500, 11) * 1e-9; % Convert to meters

% Real and Imaginary parts of the dielectric constant for silver (ε1 and ε2)
epsilon_real = [-16.0, -15.5, -14.0, -12.5, -11.5, -10.0, -9.0, -8.5, -8.0, -7.8, -7.5];
epsilon_imag = [1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0];
epsilon_complex = epsilon_real + 1i * epsilon_imag; % Combine into complex dielectric constant

% Constants
c = 3e8;  % Speed of light in vacuum (m/s)
epsilon_air = 1; % Dielectric constant of air

% Preallocate arrays for results
k_sp = zeros(1, length(wavelengths_nm));
L_sp = zeros(1, length(wavelengths_nm));
kappa_m = zeros(1, length(wavelengths_nm));
kappa_d = zeros(1, length(wavelengths_nm));

% Perform calculations for each wavelength
for i = 1:length(wavelengths_nm)
    lambda = wavelengths_nm(i); % Current wavelength in meters
    omega = 2 * pi * c / lambda; % Angular frequency
    
    % Calculate plasmon wavevector
    k_sp(i) = omega / c * sqrt((epsilon_complex(i) * epsilon_air) / (epsilon_complex(i) + epsilon_air));
    
    % Calculate propagation distance L_sp
    L_sp(i) = 1 / (2 * imag(k_sp(i)));
    
    % Calculate decay constants in metal and air
    k0 = 2 * pi / lambda;
    kappa_m(i) = k0 * sqrt(abs(epsilon_complex(i))^2 - (epsilon_complex(i) / (epsilon_complex(i) + epsilon_air)));
    kappa_d(i) = k0 * sqrt(abs(epsilon_air)^2 - (epsilon_air / (epsilon_complex(i) + epsilon_air)));
end

% Display the results
disp('Wavelength (nm)    k_sp (1/m)    L_sp (m)    kappa_m (1/m)    kappa_d (1/m)');
for i = 1:length(wavelengths_nm)
    fprintf('%9.1f        %e     %e     %e     %e\n', wavelengths_nm(i)*1e9, k_sp(i), L_sp(i), kappa_m(i), kappa_d(i));
end

% Plot results
figure;
subplot(3,1,1);
plot(wavelengths_nm*1e9, real(k_sp), 'b', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('Re(k_{sp})');
title('Real Part of Plasmon Wavevector (k_{sp})');
grid on;

subplot(3,1,2);
plot(wavelengths_nm*1e9, L_sp, 'r', 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('L_{sp} (m)');
title('Propagation Distance of Plasmons (L_{sp})');
grid on;

subplot(3,1,3);
plot(wavelengths_nm*1e9, [kappa_m; kappa_d], 'LineWidth', 2);
xlabel('Wavelength (nm)');
ylabel('\kappa (1/m)');
title('Decay Constants (\kappa_m and \kappa_d)');
legend('\kappa_m (Metal)', '\kappa_d (Air)');
grid on;

%}
% Define wavelength values in nanometers (exactly 11 values)
wavelengths_nm = linspace(500, 1500, 11) * 1e-9; % Convert to meters

% Real and Imaginary parts of the dielectric constant for silver (ε1 and ε2)
epsilon_real = [-16.0, -15.5, -14.0, -12.5, -11.5, -10.0, -9.0, -8.5, -8.0, -7.8, -7.5];
epsilon_imag = [1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0];
epsilon_complex = epsilon_real + 1i * epsilon_imag; % Combine into complex dielectric constant

% Display the dielectric constant values for silver
disp('Selected Dielectric Constant (Complex) values for Silver:');
for i = 1:length(epsilon_complex)
    fprintf('Wavelength: %4.0f nm -> Dielectric Constant: %6.2f + %.2fi\n', wavelengths_nm(i)*1e9, real(epsilon_complex(i)), imag(epsilon_complex(i)));
end

% Constants
c = 3e8;  % Speed of light in vacuum (m/s)
epsilon_air = 1; % Dielectric constant of air

% Preallocate arrays for results
k_sp = zeros(1, length(wavelengths_nm));
L_sp = zeros(1, length(wavelengths_nm));
normalized_L_sp = zeros(1, length(wavelengths_nm)); % Normalized propagation distance

% Perform calculations for each wavelength
for i = 1:length(wavelengths_nm)
    lambda = wavelengths_nm(i); % Current wavelength in meters
    omega = 2 * pi * c / lambda; % Angular frequency
    
    % Calculate plasmon wavevector
    k_sp(i) = omega / c * sqrt((epsilon_complex(i) * epsilon_air) / (epsilon_complex(i) + epsilon_air));
    
    % Calculate propagation distance L_sp
    L_sp(i) = 1 / (2 * imag(k_sp(i)));
    
    % Calculate normalized propagation distance (L_sp / lambda)
    normalized_L_sp(i) = L_sp(i) / lambda;
end

% Plot Normalized Propagation Distance vs Wavelength (Reversed)
figure;
plot(wavelengths_nm * 1e9, 1./normalized_L_sp, 'b', 'LineWidth', 2); % Invert normalized propagation distance
xlabel('Wavelength (nm)');
ylabel('Inverse of Normalized Propagation Distance (1 / L_{sp} / \lambda)');
title('Inverse of Normalized Propagation Distance of Plasmons vs Wavelength');
grid on;



%Q3 מורחבת עם הסברים

%% Question 3: Ring Resonator Analysis
%{
%% Load Data
% טוענים את הנתונים מתוך הקובץ 'data.m'
load('data_NOTUSE.m.mat');  % יש להחליף בנתיב הנכון לקובץ במחשב שלך במידת הצורך

% Assign the loaded variables to x and y
y = Loss;  % נתוני הפסדים
x = Wavelength;  % 

% נתוני הפסדים ואורכי גל, כפי שנטענו מתוך הקובץ
loss_data = Loss(8400:9100);  % הפסדים בטווח המתאים
wavelength_data = Wavelength(8400:9100);  % אורכי הגל בטווח

%% מציאת המיקום של שיא מינימום (תהודה)
[resonance_min_value, resonance_min_index] = min(loss_data);

%% חישוב עוצמת חצי הגובה
half_max_value = (max(loss_data) + resonance_min_value) / 2;

%% חישוב האינדקסים שבהם הנתונים חוצים את חצי הגובה משני צידי מינימום
left_half_index = find(loss_data(1:resonance_min_index) >= half_max_value, 1, 'last');
right_half_index = resonance_min_index + find(loss_data(resonance_min_index:end) <= half_max_value, 1, 'first') - 1;

%% חישוב רוחב ברוחב חצי גובה (FWHM)
full_width_half_max = wavelength_data(right_half_index) - wavelength_data(left_half_index);

fprintf('The FWHM is: %.2f\n', full_width_half_max);

%% ציור גרף הנתונים עם ציון נקודות FWHM
plot(wavelength_data, loss_data);
hold on;
plot(wavelength_data([left_half_index, right_half_index]), loss_data([left_half_index, right_half_index]), 'ro');
xlabel('Wavelength [μm]');
ylabel('Loss [dB]');
title('Loss Spectrum with FWHM');
legend('Data', 'FWHM Points');
hold off;

%% הגדרת lambda_resonance (אורך גל התהודה)
[~, min_index] = min(loss_data);  % מציאת המיקום של המינימום בעקומת ההפסדים
lambda_resonance = wavelength_data(min_index);  % אורך גל התהודה

%% חישוב מרווח התהודה (FSR)
[~, sorted_indices] = sort(loss_data);  % מיון לפי ערכי הפסדים
min_indices = sorted_indices(1:2);  % שני ערכי המינימום הקטנים ביותר

% חישוב מרווח התהודה בין שתי נקודות המינימום
free_spectral_range = wavelength_data(min_indices(2)) - wavelength_data(min_indices(1));

fprintf('The FSR is: %.2f\n', free_spectral_range);

%% יחס הכחדה (ER)
extinction_ratio_dB = -8.90;  % נתון ידוע בשאלה
fprintf('The ER in dB is: %.2f\n', extinction_ratio_dB);

%% חישוב Q-Factor
quality_factor = (1550) / free_spectral_range;
fprintf('The Q-factor is: %.2f\n', quality_factor);

%% חישוב n_g (מקדם השבירה הקבוצתי) – תיקון:

group_index = ((lambda_resonance*10^-9)^2) / (free_spectral_range * 10^(-9)*((2 * pi * 50 + 60*2) * 10^-6));
fprintf('The n_g is: %.2f\n', group_index);

%% חישוב alpha_dB (הפסדים ב-dB/cm) – תיקון:
L_total = (2 * pi * 50 +(2 * 60)) * 10^(-6);  % אורך המהוד הטבעתי, ודא שהיחידות נכונות
alpha_loss_dB = (4.34 * pi / (free_spectral_range * L_total)* full_width_half_max * (1 + sqrt(10^(extinction_ratio_dB / 10))));
fprintf('The alpha_dB is: %.2f\n', alpha_loss_dB*0.01);
%}

%Q4_NANO_FINAL
% Parameters
refractive_index = 2.8339;    % Effective refractive index of the material
light_wavelength = 1.55e-6;   % Light wavelength in meters
diffraction_order = -1;       % Order of diffraction (absolute value used)
entry_angle_deg = 10;         % Entry angle in degrees

% Convert the incident angle to radians
entry_angle_rad = deg2rad(entry_angle_deg);

% Grating equation: |diffraction_order| * light_wavelength = spacing * sin(entry_angle)
% Calculate the spacing of the grating
grating_spacing = abs(diffraction_order) * light_wavelength / sin(entry_angle_rad);

% Display the calculated spacing
fprintf('The calculated grating spacing is: %.6e meters\n', grating_spacing);

% Use Snell's law to determine the diffraction angle inside the material
diffraction_angle_internal_rad = asin((sin(entry_angle_rad) * refractive_index) / 1); % n_air = 1

% Convert the internal diffraction angle to degrees
diffraction_angle_deg = rad2deg(diffraction_angle_internal_rad);

% Display the diffraction angle
fprintf('The diffraction angle inside the material is: %.2f degrees\n', diffraction_angle_deg);

% Compute the grating periodicity (inverse of spacing)
grating_periodicity = 1 / grating_spacing;

% Display the grating periodicity
fprintf('The grating periodicity is: %.6e meters\n', grating_periodicity);

