%Q4

clc, clearvars

%סעיף א
% Parameters
lambda_0 = 1.55e-6;    
n_eff = 2.8339;       
n_a = 1;               
theta_a_deg = 10;      % Incident angle in degrees

% Convert angle to radians
theta_a_rad = deg2rad(theta_a_deg);

% Calculate grating periodicity (Lambda) using the provided formula
Lambda = lambda_0 / (n_eff - n_a * sin(theta_a_rad));

% Display the result
fprintf('The grating periodicity Λ is: %.6e meters\n', Lambda);

% Parameters
n_eff = 2.8339;             
n_air = 1;                  
lambda = 1.55e-6;           
m = -1;                     
theta_i_deg = 10;           

% Convert the entry angle to radians
theta_i_rad = deg2rad(theta_i_deg);

%סעיףב
% Structural data and refractive index changes
initial_n_eff = 2.8339;     % Original effective index of refraction
period = Lambda ;
delta_n_center = 0.03;      % Change at the center of the wafer
delta_n_edge = 0.06;        % Change at the edges of the wafer

% Calculate new effective indices
n_eff_center = initial_n_eff + delta_n_center;
n_eff_edge = initial_n_eff + delta_n_edge;

% calculation of the corresponding angles (subtracting lambda from n_eff)
theta_i_center = asin(n_eff_center - lambda / period);
theta_i_edge = asin(n_eff_edge - lambda / period);

% Convert radians to degrees for better visualization
theta_i_center_deg = rad2deg(theta_i_center);
theta_i_edge_deg = rad2deg(theta_i_edge);

% Display the results
disp('angle at the center (with Δn of 0.03):');
disp(theta_i_center_deg);

disp('angle at the edges (with Δn of 0.06):');
disp(theta_i_edge_deg);

%סעיף ג 
n_lambda_center= (n_eff_center-sin(theta_i_rad))*Lambda;
n_lambda_edge= (n_eff_edge-sin(theta_i_rad))*Lambda;
disp('New_lambda for center:');
disp(n_lambda_center);

disp('New_lambda for edge:');
disp(n_lambda_edge);



