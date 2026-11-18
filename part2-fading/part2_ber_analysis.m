%% BPSK BER Simulation over Rayleigh Flat Fading Channel vs. AWGN
clc; clear; close all;

% 1. Simulation Parameters
N_bits = 1e6;                  % Total number of transmitted bits
Eb = 1;                        % Bit energy (normalized to 1)
SNR_dB = 0:2:30;               % SNR (Eb/No) range in dB
SNR_linear = 10.^(SNR_dB/10);  % Convert SNR from dB to linear scale
BER_BPSK = zeros(size(SNR_dB));% Array to store simulated BER values

% Step (a): Generate random binary bit stream bk in {0, 1}
b_k = randi([0, 1], 1, N_bits);

% Step (b): BPSK Modulation (Polar NRZ: 0 -> -sqrt(Eb), 1 -> +sqrt(Eb))
x_k = (2 * b_k - 1) * sqrt(Eb);

% Loop over all SNR values (Step g)
for i = 1:length(SNR_dB)
    
    % Calculate noise power spectral density (No) for current SNR
    No = Eb / SNR_linear(i);
    
    % Step (c): Generate Rayleigh fading channel h = hr + j*hi
    % Mean = 0, Variance = 1/2 -> Std Dev = sqrt(1/2)
    h_r = sqrt(1/2) * randn(1, N_bits);
    h_i = sqrt(1/2) * randn(1, N_bits);
    h_k = h_r + 1j * h_i;
    
    % Step (c): Generate AWGN noise n = nc + j*ns
    % Mean = 0, Variance = No/2 -> Std Dev = sqrt(No/2)
    n_c = sqrt(No/2) * randn(1, N_bits);
    n_s = sqrt(No/2) * randn(1, N_bits);
    n_k = n_c + 1j * n_s;
    
    % Step (d): Compute received signal yk = hk * xk + nk
    y_k = h_k .* x_k + n_k;
    
    % Step (e): Channel Inversion (Compensate for channel gain)
    y_equalized = y_k ./ h_k;
    
    % Step (e): Hard Decision Decoding (Threshold = 0 on the real part)
    b_hat = real(y_equalized) > 0;
    
    % Step (f): Compute Bit-Error Rate (BER)
    num_errors = sum(b_k ~= b_hat);
    BER_BPSK(i) = num_errors / N_bits;
    
end

% Calculate Theoretical BER for AWGN Channel using erfc
BER_AWGN = 0.5 * erfc(sqrt(SNR_linear));

% Figure 1: BPSK over Rayleigh Flat Fading Channel Only (Step h)
figure(1);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('BER Performance of BPSK over Rayleigh Flat Fading Channel');
ylim([1e-5 1]);
legend('BPSK Rayleigh (Simulation)');

% Figure 2: Comparison between Rayleigh Fading and AWGN Channel
figure(2);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('BPSK BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]);
legend('BPSK over Rayleigh Channel', 'BPSK over AWGN Channel');

