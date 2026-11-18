%% Project Part 2: BPSK and QPSK BER Simulation over Rayleigh Flat Fading Channel
clc; clear; close all;

%% ========================================================================
%  PART 1: BPSK Simulation
%  ========================================================================

% 1. Simulation Parameters
N_bits = 1e6;                   % Total number of transmitted bits
Eb = 1;                         % Bit energy (normalized to 1)
SNR_dB = 0:2:30;                % SNR (Eb/No) range in dB
SNR_linear = 10.^(SNR_dB/10);   % Convert SNR from dB to linear scale
BER_BPSK = zeros(size(SNR_dB)); % Array to store BPSK BER values

% Step (a): Generate random binary bit stream bk in {0, 1}
b_k = randi([0, 1], 1, N_bits);

% Step (b): BPSK Modulation (0 -> -sqrt(Eb), 1 -> +sqrt(Eb))
x_BPSK = (2 * b_k - 1) * sqrt(Eb);

% Loop over all SNR values (Step g)
for i = 1:length(SNR_dB)
    No = Eb / SNR_linear(i);
    
    % Step (c): Generate Rayleigh channel (Var = 1/2) and AWGN (Var = No/2)
    h_k = sqrt(1/2) * (randn(1, N_bits) + 1j * randn(1, N_bits));
    n_k = sqrt(No/2) * (randn(1, N_bits) + 1j * randn(1, N_bits));
    
    % Step (d): Received signal
    y_k = h_k .* x_BPSK + n_k;
    
    % Step (e): Channel Inversion and Hard Decision (Threshold = 0)
    y_eq = y_k ./ h_k;
    b_hat_BPSK = real(y_eq) > 0;
    
    % Step (f): Compute BPSK BER
    BER_BPSK(i) = sum(b_k ~= b_hat_BPSK) / N_bits;
end

% Theoretical BER for BPSK/QPSK over AWGN Channel
BER_AWGN = 0.5 * erfc(sqrt(SNR_linear));

% Figure 1: BPSK over Rayleigh Channel Only
figure(1);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('BER Performance of BPSK over Rayleigh Flat Fading Channel');
ylim([1e-5 1]);
legend('BPSK Rayleigh (Simulation)');

% Figure 2: BPSK Comparison (Rayleigh vs. AWGN)
figure(2);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('BPSK BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]);
legend('BPSK over Rayleigh Channel', 'BPSK over AWGN Channel');


%% ========================================================================
%  PART 2: QPSK Simulation
%  ========================================================================

N_symbols_QPSK = N_bits / 2;    % Each QPSK symbol carries 2 bits
BER_QPSK = zeros(size(SNR_dB)); % Array to store QPSK BER values

% Step (b): QPSK Modulation based on given constellation (b0 -> I, b1 -> Q)
b_0 = b_k(1:2:end);             % First bit of each symbol (odd indices)
b_1 = b_k(2:2:end);             % Second bit of each symbol (even indices)

I_QPSK = (2 * b_0 - 1) * sqrt(Eb);
Q_QPSK = (2 * b_1 - 1) * sqrt(Eb);
x_QPSK = I_QPSK + 1j * Q_QPSK;  % Complex QPSK symbols (Es = 2*Eb)

% Loop over all SNR values
for i = 1:length(SNR_dB)
    No = Eb / SNR_linear(i);
    
    % Step (c): Generate Rayleigh channel and AWGN of length N_symbols_QPSK
    h_k = sqrt(1/2) * (randn(1, N_symbols_QPSK) + 1j * randn(1, N_symbols_QPSK));
    n_k = sqrt(No/2) * (randn(1, N_symbols_QPSK) + 1j * randn(1, N_symbols_QPSK));
    
    % Step (d): Received signal
    y_k = h_k .* x_QPSK + n_k;
    
    % Step (e): Channel Inversion and Hard Decision on I and Q (Threshold = 0)
    y_eq = y_k ./ h_k;
    b0_hat = real(y_eq) > 0;
    b1_hat = imag(y_eq) > 0;
    
    % Reconstruct estimated bit stream
    b_hat_QPSK = zeros(1, N_bits);
    b_hat_QPSK(1:2:end) = b0_hat;
    b_hat_QPSK(2:2:end) = b1_hat;
    
    % Step (f): Compute QPSK BER
    BER_QPSK(i) = sum(b_k ~= b_hat_QPSK) / N_bits;
end

% Figure 3: QPSK over Rayleigh Channel Only
figure(3);
semilogy(SNR_dB, BER_QPSK, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('BER Performance of QPSK over Rayleigh Flat Fading Channel');
ylim([1e-5 1]);
legend('QPSK Rayleigh (Simulation)');

% Figure 4: QPSK Comparison (Rayleigh vs. AWGN)
figure(4);
semilogy(SNR_dB, BER_QPSK, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on;
xlabel('SNR (E_b/N_0) [dB]');
ylabel('Bit Error Rate (BER)');
title('QPSK BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]);
legend('QPSK over Rayleigh Channel', 'QPSK over AWGN Channel');

