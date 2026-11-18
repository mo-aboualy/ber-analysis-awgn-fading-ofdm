%% Digital Communications Project Part 2
%  BER Performance for BPSK, QPSK, and 16-QAM over Rayleigh Flat Fading Channel
clc; clear; close all;

% Common Simulation Parameters
N_bits = 1e6;                   % Total number of bits (divisible by 4 for 16-QAM)
SNR_dB = 0:2:30;                % SNR (Eb/No) range in dB
SNR_linear = 10.^(SNR_dB/10);   % Linear SNR values

% Step (a): Generate random binary bit stream bk in {0, 1}
b_k = randi([0, 1], 1, N_bits);

%% ========================================================================
%  PART 1: BPSK Simulation
%  ========================================================================
Eb_BPSK = 1;
BER_BPSK = zeros(size(SNR_dB));

% Step (b): BPSK Modulation (0 -> -1, 1 -> +1)
x_BPSK = (2 * b_k - 1) * sqrt(Eb_BPSK);

for i = 1:length(SNR_dB)
    No = Eb_BPSK / SNR_linear(i);
    
    % Steps (c, d, e): Pass through Rayleigh channel and apply channel inversion
    y_eq = pass_through_rayleigh_channel(x_BPSK, No);
    
    % Step (e): Hard Decision Decoding (Threshold = 0)
    b_hat_BPSK = real(y_eq) > 0;
    
    % Step (f): Compute BPSK BER
    BER_BPSK(i) = sum(b_k ~= b_hat_BPSK) / N_bits;
end

% Theoretical BER for BPSK/QPSK over AWGN Channel
BER_AWGN_BPSK_QPSK = 0.5 * erfc(sqrt(SNR_linear));

% Figure 1: BPSK over Rayleigh Channel Only
figure(1);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('BER Performance of BPSK over Rayleigh Flat Fading Channel');
ylim([1e-5 1]); legend('BPSK Rayleigh (Simulation)');

% Figure 2: BPSK Comparison (Rayleigh vs. AWGN)
figure(2);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN_BPSK_QPSK, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('BPSK BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]); legend('BPSK over Rayleigh Channel', 'BPSK over AWGN Channel');

%% ========================================================================
%  PART 2: QPSK Simulation
%  ========================================================================
Eb_QPSK = 1;                    % Es = 1^2 + 1^2 = 2 -> Eb = Es/2 = 1
BER_QPSK = zeros(size(SNR_dB));

% Step (b): QPSK Modulation (b0 -> I, b1 -> Q)
b_0_qpsk = b_k(1:2:end);
b_1_qpsk = b_k(2:2:end);
x_QPSK = (2 * b_0_qpsk - 1) + 1j * (2 * b_1_qpsk - 1);

for i = 1:length(SNR_dB)
    No = Eb_QPSK / SNR_linear(i);
    
    % Steps (c, d, e): Pass through Rayleigh channel and apply channel inversion
    y_eq = pass_through_rayleigh_channel(x_QPSK, No);
    
    % Step (e): Hard Decision Decoding on I and Q (Threshold = 0)
    b_hat_QPSK = zeros(1, N_bits);
    b_hat_QPSK(1:2:end) = real(y_eq) > 0;
    b_hat_QPSK(2:2:end) = imag(y_eq) > 0;
    
    % Step (f): Compute QPSK BER
    BER_QPSK(i) = sum(b_k ~= b_hat_QPSK) / N_bits;
end

% Figure 3: QPSK over Rayleigh Channel Only
figure(3);
semilogy(SNR_dB, BER_QPSK, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('BER Performance of QPSK over Rayleigh Flat Fading Channel');
ylim([1e-5 1]); legend('QPSK Rayleigh (Simulation)');

% Figure 4: QPSK Comparison (Rayleigh vs. AWGN)
figure(4);
semilogy(SNR_dB, BER_QPSK, 'm-d', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN_BPSK_QPSK, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('QPSK BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]); legend('QPSK over Rayleigh Channel', 'QPSK over AWGN Channel');

%% ========================================================================
%  PART 3: 16-QAM Simulation
%  ========================================================================
Eb_16QAM = 2.5;                 % Es = 10 -> Eb = Es/4 = 2.5
BER_16QAM = zeros(size(SNR_dB));

% Extract 4 bits per symbol (b0, b1, b2, b3)
b_0 = b_k(1:4:end);
b_1 = b_k(2:4:end);
b_2 = b_k(3:4:end);
b_3 = b_k(4:4:end);

% Step (b): 16-QAM Gray Mapping based on given constellation
% 00 -> -3, 01 -> -1, 11 -> +1, 10 -> +3
I_16QAM = (2 * b_0 - 1) .* (3 - 2 * b_1);
Q_16QAM = (2 * b_2 - 1) .* (3 - 2 * b_3);
x_16QAM = I_16QAM + 1j * Q_16QAM;

for i = 1:length(SNR_dB)
    No = Eb_16QAM / SNR_linear(i);
    
    % Steps (c, d, e): Pass through Rayleigh channel and apply channel inversion
    y_eq = pass_through_rayleigh_channel(x_16QAM, No);
    
    % Step (e): Hard Decision Decoding for 16-QAM (Thresholds at 0 and +/-2)
    I_rx = real(y_eq);
    Q_rx = imag(y_eq);
    
    b0_hat = I_rx > 0;          % Sign of I determines b0
    b1_hat = abs(I_rx) < 2;     % Magnitude of I (<2 means +/-1, >2 means +/-3)
    b2_hat = Q_rx > 0;          % Sign of Q determines b2
    b3_hat = abs(Q_rx) < 2;     % Magnitude of Q (<2 means +/-1, >2 means +/-3)
    
    % Reconstruct estimated bit stream
    b_hat_16QAM = zeros(1, N_bits);
    b_hat_16QAM(1:4:end) = b0_hat;
    b_hat_16QAM(2:4:end) = b1_hat;
    b_hat_16QAM(3:4:end) = b2_hat;
    b_hat_16QAM(4:4:end) = b3_hat;
    
    % Step (f): Compute 16-QAM BER
    BER_16QAM(i) = sum(b_k ~= b_hat_16QAM) / N_bits;
end

% Theoretical BER for 16-QAM over AWGN Channel
BER_AWGN_16QAM = (3/8) * erfc(sqrt((2/5) * SNR_linear));

% Figure 5: 16-QAM over Rayleigh Channel Only
figure(5);
semilogy(SNR_dB, BER_16QAM, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('BER Performance of 16-QAM over Rayleigh Flat Fading Channel');
ylim([1e-5 1]); legend('16-QAM Rayleigh (Simulation)');

% Figure 6: 16-QAM Comparison (Rayleigh vs. AWGN)
figure(6);
semilogy(SNR_dB, BER_16QAM, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_AWGN_16QAM, 'r--s', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('16-QAM BER Comparison: Rayleigh Flat Fading vs. AWGN Channel');
ylim([1e-5 1]); legend('16-QAM over Rayleigh Channel', '16-QAM over AWGN Channel');

%% ========================================================================
%  Overall Comparison Figure: BPSK vs. QPSK vs. 16-QAM over Rayleigh
%  ========================================================================
figure(7);
semilogy(SNR_dB, BER_BPSK, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
semilogy(SNR_dB, BER_QPSK, 'm--d', 'LineWidth', 1.5, 'MarkerSize', 6);
semilogy(SNR_dB, BER_16QAM, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6);
grid on; xlabel('SNR (E_b/N_0) [dB]'); ylabel('Bit Error Rate (BER)');
title('BER Comparison: BPSK, QPSK, and 16-QAM over Rayleigh Channel');
ylim([1e-5 1]);
legend('BPSK', 'QPSK', '16-QAM');

%% ========================================================================
%  Helper Function: Rayleigh Channel, AWGN, and Channel Inversion
%  ========================================================================
function y_eq = pass_through_rayleigh_channel(x_k, No)
    N_sym = length(x_k);
    
    % Step (c): Generate Rayleigh fading channel (Mean = 0, Var = 1/2)
    h_k = sqrt(1/2) * (randn(1, N_sym) + 1j * randn(1, N_sym));
    
    % Step (c): Generate AWGN noise (Mean = 0, Var = No/2)
    n_k = sqrt(No/2) * (randn(1, N_sym) + 1j * randn(1, N_sym));
    
    % Step (d): Compute received symbol vector yk = hk * xk + nk
    y_k = h_k .* x_k + n_k;
    
    % Step (e): Channel Inversion (compensate for channel gain)
    y_eq = y_k ./ h_k;
end

