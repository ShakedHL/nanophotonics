%Q2 plasmon

clc, clearvars

wavelengths_nm = linspace(500, 1500, 11) * 1e-9; 

epsilon_real = [-16.0, -15.5, -14.0, -12.5, -11.5, -10.0, -9.0, -8.5, -8.0, -7.8, -7.5];
epsilon_imag = [1.5, 1.8, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0];
epsilon_complex = epsilon_real + 1i * epsilon_imag; 

% Display the dielectric constant values for silver
disp('Selected Dielectric Constant (Complex) values for Silver:');
for i = 1:length(epsilon_complex)
    fprintf('Wavelength: %4.0f nm -> Dielectric Constant: %6.2f + %.2fi\n', wavelengths_nm(i)*1e9, real(epsilon_complex(i)), imag(epsilon_complex(i)));
end

c = 3e8;  
epsilon_air = 1; 

k_sp = zeros(1, length(wavelengths_nm));
L_sp = zeros(1, length(wavelengths_nm));
normalized_L_sp = zeros(1, length(wavelengths_nm)); 

for i = 1:length(wavelengths_nm)
    lambda = wavelengths_nm(i); 
    omega = 2 * pi * c / lambda; 
    
    k_sp(i) = omega / c * sqrt((epsilon_complex(i) * epsilon_air) / (epsilon_complex(i) + epsilon_air));
    L_sp(i) = 1 / (2 * imag(k_sp(i)));
    normalized_L_sp(i) = L_sp(i) / lambda;
end

figure;
plot(wavelengths_nm * 1e9, 1./normalized_L_sp, 'b', 'LineWidth', 2); 
xlabel('Wavelength (nm)');
ylabel('Inverse of Normalized Propagation Distance (1 / L_{sp} / \lambda)');
title('Inverse of Normalized Propagation Distance of Plasmons vs Wavelength');
grid on;


