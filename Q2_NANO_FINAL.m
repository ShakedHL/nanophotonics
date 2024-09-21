%Q2 

% Clear workspace and close all figures
clear;
close all;


% Define the wavelength range (in nm)
selected_wavelengths = 500:50:1000;

% Dielectric constants for silver at the selected wavelengths (data from Palik: https://refractiveindex.info/?shelf=main&book=Ag&page=Ferrera-298K)
epsilon_silver = [
    -9.9801 + 0.26143i;  % 500 nm
    -13.031 + 0.31637i;  % 550 nm
    -16.335 + 0.38326i;  % 600 nm
    -19.891 + 0.45417i;  % 650 nm
    -23.705 + 0.53424i;  % 700 nm
    -27.816 + 0.63413i;  % 750 nm
    -32.237 + 0.72624i;  % 800 nm
    -36.919 + 0.85236i;  % 850 nm
    -41.845 + 0.96771i;  % 900 nm
    -47.123 + 1.1135i;   % 950 nm
    -52.615 + 1.2771i;   % 1000 nm
];

% Air dielectric constant
epsilon_air = 1;

% Preallocate arrays for decay constants and penetration depths
decay_silver = zeros(1, length(selected_wavelengths));
decay_air = zeros(1, length(selected_wavelengths));
plasmon_wavelength_air = zeros(1, length(selected_wavelengths));
plasmon_wavelength_silver = zeros(1, length(selected_wavelengths));
L_spp = zeros(1, length(selected_wavelengths));
beta2 = zeros(1, length(selected_wavelengths));
penetration_depth_air = zeros(1, length(selected_wavelengths));
penetration_depth_silver = zeros(1, length(selected_wavelengths));
plasmon_wavelength = zeros(1, length(selected_wavelengths));

% Loop over each wavelength to compute the relevant parameters
for idx = 1:length(selected_wavelengths)
    lambda = selected_wavelengths(idx) * 1e-9; % Convert to meters
    epsilon_metal = epsilon_silver(idx);
    k0 = 2 * pi / lambda; % Free space wavenumber
    beta(idx) = k0 * sqrt((epsilon_air * real(epsilon_metal)) / (epsilon_air + real(epsilon_metal)));
    beta2(idx) = (k0 * epsilon_air * imag(epsilon_metal) / (2 * real(epsilon_metal)^2)) * (beta(idx) / k0)^3;
    
    % Calculate the surface plasmon propagation length
    L_spp(idx) = 1 / (2 * beta2(idx));
    
    % Calculate the penetration depths in air and silver
    penetration_depth_air(idx) = abs((1 / k0) * sqrt(abs(real(epsilon_metal) + epsilon_air) / (epsilon_air)^2));
    penetration_depth_silver(idx) = abs((1 / k0) * sqrt(abs(real(epsilon_metal) + epsilon_air) / (epsilon_metal)^2));
    
    % Plasmon wavelength calculation
    plasmon_wavelength(idx) = (2 * pi) / real(beta(idx));
end

% Plot of normalized progression distance vs wavelength
if any(L_spp)
    figure(1);
    plot(selected_wavelengths, L_spp);
    xlabel('Wavelength (nm)');
    ylabel('Progression Distance');
    title('Progression Distance vs Wavelength');
    grid on;
end

% Plot of decay constant vs wavelength
if any(beta2)
    figure(2);
    plot(selected_wavelengths, beta2);
    xlabel('Wavelength (nm)');
    ylabel('Decay Constant');
    title('Decay Constant vs Wavelength');
    grid on;
end

% Plot penetration depth in air and silver
if any(penetration_depth_air) && any(penetration_depth_silver)
    figure(3);
    subplot(2, 1, 1);
    plot(selected_wavelengths, penetration_depth_air);
    xlabel('Wavelength (nm)');
    ylabel('Penetration Depth (m)');
    title('Penetration Depth in Air vs Wavelength');
    legend('Air');
    grid on;

    subplot(2, 1, 2);
    plot(selected_wavelengths, penetration_depth_silver);
    xlabel('Wavelength (nm)');
    ylabel('Penetration Depth (m)');
    title('Penetration Depth in Silver vs Wavelength');
    legend('Silver');
    grid on;
end

% Plot normalized progression distance vs wavelength
if any(L_spp)
    figure(4);
    plot(selected_wavelengths, L_spp / lambda);
    xlabel('Wavelength (nm)');
    ylabel('Normalized Progression Distance');
    title('Normalized Progression Distance vs Wavelength');
    grid on;
end

% Create table of results and display
results = table(selected_wavelengths', plasmon_wavelength' * 1e9, beta' * 1e6, beta2' * 1e6, L_spp' * 1e6, penetration_depth_air' * 1e6, penetration_depth_silver' * 1e6, ...
    'VariableNames', {'Wavelength (nm)', 'Plasmon Wavelength (nm)', 'Beta (micron)', 'Decay Constant (micron)', 'L_spp (micron)', 'Penetration Depth Air (micron)', 'Penetration Depth Silver (micron)'});

disp('Calculation Results:');
disp(results);
