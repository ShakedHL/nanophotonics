clc, clearvars

%% Question 3: Ring Resonator 

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

