clc, clearvars

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
