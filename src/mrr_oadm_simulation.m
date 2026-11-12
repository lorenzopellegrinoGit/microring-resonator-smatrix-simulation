%% POC Assignment 5 - Lorenzo Pellegrino
clc;
clearvars;
close all;

%% Parametri globali
neff0 = 2.6;                                      % Effective refractive index
lambda_span = linspace(1.545e-6, 1.565e-6, 2000); % Wavelength vector [m]
alpha_loss = 0;                                   % Ring loss coefficient [1/m]
S41_mag = 0.05;                                   % Non-ideality reflection parameter

%% --- SEZIONE Single MRR
r1 = 8e-6;                                        % Base ring radius [m]
r2 = 10e-6;                                       % Larger ring radius [m]

k_pow1 = 0.10;                                    % Base power coupling coefficient
k_pow2 = 0.25;                                    % Higher power coupling coefficient

[T_pow1, R_pow1] = simulate_single_MRR(lambda_span, r1, neff0, k_pow1, S41_mag, alpha_loss);
[T_pow2, R_pow2] = simulate_single_MRR(lambda_span, r1, neff0, k_pow2, S41_mag, alpha_loss);
[T_pow3, R_pow3] = simulate_single_MRR(lambda_span, r2, neff0, k_pow1, S41_mag, alpha_loss);

figure(1);
set(gcf, 'Position', [100, 100, 1000, 450]);

subplot(1,2,1);
plot(lambda_span * 1e6, T_pow1, 'b', 'LineWidth', 1.5); hold on;
plot(lambda_span * 1e6, T_pow2, 'g--', 'LineWidth', 1.5);
plot(lambda_span * 1e6, T_pow3, 'r:', 'LineWidth', 1.5);
xlabel('Wavelength (\mum)');
ylabel('|T|^2 (Normalized Power)');
title('Transmission Function (Through Port)');
legend(sprintf('r = 8\\mum, K = %.2f', k_pow1), ...
       sprintf('r = 8\\mum, K = %.2f', k_pow2), ...
       sprintf('r = 10\\mum, K = %.2f', k_pow1), 'Location', 'SouthEast');
grid on;

subplot(1,2,2);
plot(lambda_span * 1e6, R_pow1, 'b', 'LineWidth', 1.5); hold on;
plot(lambda_span * 1e6, R_pow2, 'g--', 'LineWidth', 1.5);
plot(lambda_span * 1e6, R_pow3, 'r:', 'LineWidth', 1.5);
xlabel('Wavelength (\mum)');
ylabel('|\Gamma_{in}|^2 (Normalized Power)');
title('Input Reflectivity (Return Loss)');
legend(sprintf('r = 8\\mum, K = %.2f', k_pow1), ...
       sprintf('r = 8\\mum, K = %.2f', k_pow2), ...
       sprintf('r = 10\\mum, K = %.2f', k_pow1), 'Location', 'NorthEast');
grid on;

sgtitle(sprintf('Single Micro-Ring Resonator (MRR) - S_{41} = %.2f', S41_mag));

%% --- SEZIONE Optical Add-Drop Multiplexer (OADM)
k1_pow_A = 0.10; k2_pow_A = 0.10;                 % Power coupling coefficients (Design A)
k1_pow_B = 0.20; k2_pow_B = 0.20;                 % Power coupling coefficients (Design B)

[P_transA, P_dropA] = simulate_OADM(lambda_span, r1, neff0, k1_pow_A, k2_pow_A, alpha_loss);
[P_transB, P_dropB] = simulate_OADM(lambda_span, r1, neff0, k1_pow_B, k2_pow_B, alpha_loss);

figure(2);
set(gcf, 'Position', [150, 150, 1000, 450]);

subplot(1,2,1);
plot(lambda_span * 1e6, P_transA, 'g', 'LineWidth', 1.5); hold on;
plot(lambda_span * 1e6, P_transB, 'b--', 'LineWidth', 1.5);
xlabel('Wavelength (\mum)');
ylabel('I_{trans} / I_{in}');
title('Transmission (Through Port)');
legend('K_1=K_2=0.10', 'K_1=K_2=0.20', 'Location', 'SouthEast');
ylim([0 1.05]);
grid on;

subplot(1,2,2);
plot(lambda_span * 1e6, P_dropA, 'm', 'LineWidth', 1.5); hold on;
plot(lambda_span * 1e6, P_dropB, 'r--', 'LineWidth', 1.5);
xlabel('Wavelength (\mum)');
ylabel('I_{drop} / I_{in}');
title('Drop Function (Drop Port)');
legend('K_1=K_2=0.10', 'K_1=K_2=0.20', 'Location', 'NorthEast');
ylim([0 1.05]);
grid on;

sgtitle('OADM Ring Resonator - Different Coupling Coefficients');

%% --- FSR and Q FACTOR
L_base = 2 * pi * r1;                             % Ring circumference [m]
a_base = exp(-alpha_loss * L_base / 2);           % Round-trip amplitude attenuation
t1_base = sqrt(1 - k1_pow_A);                     % Direct transmission coefficient 1
t2_base = sqrt(1 - k2_pow_A);                     % Direct transmission coefficient 2

[~, idx_res] = max(P_dropA);
lambda_res = lambda_span(idx_res);                % Resonance wavelength [m]

FSR_teo = (lambda_res^2) / (neff0 * L_base);      % Theoretical FSR [m]
FWHM_teo = ((1 - t1_base * t2_base * a_base) * lambda_res^2) / (pi * neff0 * L_base * sqrt(t1_base * t2_base * a_base)); % Theoretical FWHM [m]
Q_teo = lambda_res / FWHM_teo;                    % Quality factor [dimensionless]

fprintf('\n--- OADM PARAMETER ANALYSIS RESULTS (r=8um, K1=K2=0.10) ---\n');
fprintf('Resonance Wavelength         : %.3f nm\n', lambda_res * 1e9);
fprintf('Free Spectral Range (FSR)    : %.2f nm\n', FSR_teo * 1e9);
fprintf('Full Width at Half Max (FWHM): %.2f pm\n', FWHM_teo * 1e12);
fprintf('Quality Factor (Q)           : %.0f\n', Q_teo);
fprintf('-----------------------------------------------------------\n\n');

%% --- Funzioni locali
function [T_power, R_power] = simulate_single_MRR(lambda_span, r, neff, k_power, S41_mag, alpha_loss)
    L = 2 * pi * r;                               % Ring circumference [m]
    a = exp(-alpha_loss * L / 2);                 % Amplitude attenuation
    C = sqrt(k_power);                            % Cross-coupling amplitude
    t = sqrt(1 - C^2 - S41_mag^2);                % Direct transmission amplitude
    
    S21 = t;
    S31 = 1i * C;           
    S42 = 1i * C;           
    S41 = S41_mag;          
    
    theta = (2 * pi * neff * L) ./ lambda_span;   % Phase shift per round-trip [rad]
    feedback_term = a .* exp(-1i * theta);        % Round-trip feedback term
    
    Gamma_in = (2 * S41 * S42 .* feedback_term) ./ (1 - S21 .* feedback_term);
    T = S21 + ((S31^2 + S41^2) .* feedback_term) ./ (1 - S21 .* feedback_term);
    
    R_power = abs(Gamma_in).^2;                   % Normalized reflected power
    T_power = abs(T).^2;                          % Normalized transmitted power
end

function [P_trans, P_drop] = simulate_OADM(lambda_span, r, neff, k1_pow, k2_pow, alpha_loss)
    L = 2 * pi * r;                               % Ring circumference [m]
    a = exp(-alpha_loss * L / 2);                 % Amplitude attenuation
    theta = (2 * pi * neff * L) ./ lambda_span;   % Phase shift per round-trip [rad]
    
    t1 = sqrt(1 - k1_pow);                        % Direct transmission amplitude (Coupler 1)
    C1 = sqrt(k1_pow);                            % Cross-coupling amplitude (Coupler 1)
    t2 = sqrt(1 - k2_pow);                        % Direct transmission amplitude (Coupler 2)
    C2 = sqrt(k2_pow);                            % Cross-coupling amplitude (Coupler 2)
    
    feedback_OADM = a .* exp(-1i * theta);        % Round-trip feedback term
    
    E_trans = (t1 - t2 .* feedback_OADM) ./ (1 - t1 .* t2 .* feedback_OADM);
    E_drop  = (-C1 * C2 .* sqrt(a) .* exp(-1i * theta / 2)) ./ (1 - t1 .* t2 .* feedback_OADM);
    
    P_trans = abs(E_trans).^2;                    % Normalized transmitted power
    P_drop  = abs(E_drop).^2;                     % Normalized dropped power
end