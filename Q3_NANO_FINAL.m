clc, clearvars

%% Question 3: Ring Resonator Analysis

%% Load Data
load('data_NOTUSE.m.mat');

% Assign the loaded variables to x and y
y = Loss;
x = Wavelength;

% Define the range of data
loss_data = Loss(8400:9100);
wavelength_data = Wavelength(8400:9100);

%% Find the minimum peak (Resonance)
[resonance_min_value, resonance_min_index] = min(loss_data);

half_max_value = (max(loss_data) + resonance_min_value) / 2;

%% Find the indices that cross the half-maximum value
left_half_index = find(loss_data(1:resonance_min_index) >= half_max_value, 1, 'last');
right_half_index = resonance_min_index + find(loss_data(resonance_min_index:end) <= half_max_value, 1, 'first') - 1;

%% Calculate Full Width at Half Maximum (FWHM)
full_width_half_max = wavelength_data(right_half_index) - wavelength_data(left_half_index);

fprintf('The FWHM is: %.2f\n', full_width_half_max);

%% Plot the data and FWHM points
plot(wavelength_data, loss_data);
hold on;
plot(wavelength_data([left_half_index, right_half_index]), loss_data([left_half_index, right_half_index]), 'ro');
xlabel('Wavelength [μm]');
ylabel('Loss [dB]');
title('Loss Spectrum with FWHM');
legend('Data', 'FWHM Points');
hold off;

%% Calculate lambda_resonance
[~, min_index] = min(loss_data);
lambda_resonance = wavelength_data(min_index);

%% Calculate Free Spectral Range (FSR)
[~, sorted_indices] = sort(loss_data);
min_indices = sorted_indices(1:2);
free_spectral_range = wavelength_data(min_indices(2)) - wavelength_data(min_indices(1));

fprintf('The FSR is: %.2f\n', free_spectral_range);

%% Extinction Ratio (ER)
extinction_ratio_dB = -8.90;
fprintf('The ER in dB is: %.2f\n', extinction_ratio_dB);

%% Calculate Q-Factor
quality_factor = (1550) / free_spectral_range;
fprintf('The Q-factor is: %.2f\n', quality_factor);

%% Calculate Group Index (n_g)
group_index = ((lambda_resonance * 10^-9)^2) / (free_spectral_range * 10^(-9) * ((2 * pi * 50 + 60 * 2) * 10^-6));
fprintf('The n_g is: %.2f\n', group_index);

%% Calculate Loss in dB/cm (alpha_dB)
L_total = (2 * pi * 50 + (2 * 60)) * 10^(-6);
alpha_loss_dB = (4.34 * pi / (free_spectral_range * L_total) * full_width_half_max * (1 + sqrt(10^(extinction_ratio_dB / 10))));
fprintf('The alpha_dB is: %.2f\n', alpha_loss_dB * 0.01);
