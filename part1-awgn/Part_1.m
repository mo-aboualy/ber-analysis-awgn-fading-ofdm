%% 0. Environment Initialization
clearvars;
clc;
close all;
rng('default'); % Ensure reproducibility of Monte Carlo simulation results


%% ========================================================================
%  1. SIMULATION PARAMETERS & INITIALIZATION
%  ========================================================================

% Number of transmitted bits (chosen as a multiple of 12 to be evenly 
% divisible by log2(M) = 1, 2, 3, and 4 bits/symbol without padding)
numBits             = 1.2e6;

% Signal-to-Noise Ratio (Eb/No) parameters
EbNo_dB             = -2 : 1 : 14;
EbNo_lin            = 10 .^ (EbNo_dB / 10);
numSNR              = length(EbNo_dB);

% Generate uniformly distributed random binary data stream
txBits              = randi([0, 1], 1, numBits);

% Pre-allocate memory for simulated BER vectors
BER_sim_BPSK        = zeros(1, numSNR);
BER_sim_QPSK_Gray   = zeros(1, numSNR);
BER_sim_QPSK_NoGray = zeros(1, numSNR);
BER_sim_8PSK        = zeros(1, numSNR);
BER_sim_16QAM       = zeros(1, numSNR);
BER_sim_BFSK        = zeros(1, numSNR);


%% ========================================================================
%  2. BASEBAND MAPPER (TRANSMITTER)
%  ========================================================================

% -------------------------------------------------------------------------
% 2.1 BPSK Mapper (Figure 2)
%     Mapping : '0' -> -1, '1' -> +1
%     Energy  : Es = Eb = 1
% -------------------------------------------------------------------------
txSym_BPSK          = 2 * txBits - 1;

% -------------------------------------------------------------------------
% 2.2 QPSK Mapper - Gray Coded (Figure 2: b0 b1)
%     Mapping : 00 -> -1-1j | 01 -> -1+1j | 11 -> +1+1j | 10 -> +1-1j
%     Energy  : Es = 2  =>  Eb = Es / 2 = 1
% -------------------------------------------------------------------------
bits_QPSK           = reshape(txBits, 2, []).';          % [N/2 x 2] matrix
I_QPSK_Gray         = 2 * bits_QPSK(:, 1) - 1;           % In-phase (b0)
Q_QPSK_Gray         = 2 * bits_QPSK(:, 2) - 1;           % Quadrature (b1)
txSym_QPSK_Gray     = I_QPSK_Gray + 1j * Q_QPSK_Gray;

% -------------------------------------------------------------------------
% 2.3 QPSK Mapper - Non-Gray Coded (Figure 3: b0 b1)
%     Mapping : 00 -> -1-1j | 01 -> -1+1j | 10 -> +1+1j | 11 -> +1-1j
%     Logic   : I-component is determined by b0
%               Q-component is determined by (b0 XOR b1)
%     Energy  : Es = 2  =>  Eb = Es / 2 = 1
% -------------------------------------------------------------------------
I_QPSK_NoGray       = 2 * bits_QPSK(:, 1) - 1;
Q_QPSK_NoGray       = 2 * xor(bits_QPSK(:, 1), bits_QPSK(:, 2)) - 1;
txSym_QPSK_NoGray   = I_QPSK_NoGray + 1j * Q_QPSK_NoGray;

% -------------------------------------------------------------------------
% 2.4 8PSK Mapper - Gray Coded (Figure 2: b0 b1 b2)
%     Mapping : Unit-circle phases indexed by decimal value of [b0 b1 b2]
%               000 (0) ->   0     | 001 (1) ->  pi/4
%               010 (2) -> 3*pi/4  | 011 (3) ->  pi/2
%               100 (4) -> -pi/4   | 101 (5) -> -pi/2
%               110 (6) ->  pi     | 111 (7) -> -3*pi/4
%     Energy  : Es = 1  =>  Eb = Es / 3 = 1/3
% -------------------------------------------------------------------------
phaseAngles_8PSK    = [0, pi/4, 3*pi/4, pi/2, -pi/4, -pi/2, pi, -3*pi/4];
constellation_8PSK  = exp(1j * phaseAngles_8PSK).';      % [8 x 1] LUT

bits_8PSK           = reshape(txBits, 3, []).';          % [N/3 x 3] matrix
symIdx_8PSK         = bits_8PSK(:, 1) * 4 + ...
                      bits_8PSK(:, 2) * 2 + ...
                      bits_8PSK(:, 3) + 1;               % 1-based index
txSym_8PSK          = constellation_8PSK(symIdx_8PSK);

% -------------------------------------------------------------------------
% 2.5 16-QAM Mapper - Gray Coded (Figure 2: b0 b1 b2 b3)
%     Mapping : I-axis (b0 b1) & Q-axis (b2 b3)
%               00 -> -3 | 01 -> -1 | 11 -> +1 | 10 -> +3
%     Energy  : Es = 10  =>  Eb = Es / 4 = 2.5
% -------------------------------------------------------------------------
bits_16QAM          = reshape(txBits, 4, []).';          % [N/4 x 4] matrix
I_16QAM             = (2 * bits_16QAM(:, 1) - 1) .* (3 - 2 * bits_16QAM(:, 2));
Q_16QAM             = (2 * bits_16QAM(:, 3) - 1) .* (3 - 2 * bits_16QAM(:, 4));
txSym_16QAM         = I_16QAM + 1j * Q_16QAM;

% -------------------------------------------------------------------------
% 2.6 BFSK Mapper - Orthogonal Signal-Space Representation (Task 1.5)
%     Basis   : phi_1(t) = sqrt(2/Tb)*cos(2*pi*f1*t)
%               phi_2(t) = sqrt(2/Tb)*cos(2*pi*f2*t)
%     Mapping : Bit '0' (s1) -> [sqrt(Eb), 0]  (1 + 0j)
%               Bit '1' (s2) -> [0, sqrt(Eb)]  (0 + 1j)
%     Energy  : Es = Eb = 1
% -------------------------------------------------------------------------
Eb_BFSK             = 1.0;
phi1_BFSK           = (txBits == 0) * sqrt(Eb_BFSK);     % Projection on phi_1
phi2_BFSK           = (txBits == 1) * sqrt(Eb_BFSK);     % Projection on phi_2
txSym_BFSK          = phi1_BFSK + 1j * phi2_BFSK;


%% ========================================================================
%  3. AWGN CHANNEL & MAXIMUM LIKELIHOOD DEMAPPER LOOP
%  ========================================================================

fprintf('==============================================================================\n');
fprintf('  Starting Baseband Monte Carlo BER Simulation (%d bits)\n', numBits);
fprintf('==============================================================================\n');

for idx = 1 : numSNR
    
    current_EbNo = EbNo_lin(idx);
    
    % ---------------------------------------------------------------------
    % 3.1 BPSK Channel & Demapper (Eb = 1)
    % ---------------------------------------------------------------------
    Eb_BPSK             = 1.0;
    No_BPSK             = Eb_BPSK / current_EbNo;
    noise_BPSK          = sqrt(No_BPSK / 2) * randn(size(txSym_BPSK));
    
    rxSym_BPSK          = txSym_BPSK + noise_BPSK;
    rxBits_BPSK         = (rxSym_BPSK > 0);
    
    BER_sim_BPSK(idx)   = sum(rxBits_BPSK ~= txBits) / numBits;
    
    % ---------------------------------------------------------------------
    % 3.2 QPSK Channel & Demapper - Gray & Non-Gray (Eb = 1)
    % ---------------------------------------------------------------------
    Eb_QPSK             = 1.0;
    No_QPSK             = Eb_QPSK / current_EbNo;
    noise_QPSK          = sqrt(No_QPSK / 2) * ...
                          (randn(size(txSym_QPSK_Gray)) + 1j * randn(size(txSym_QPSK_Gray)));
    
    % --- Case A: Figure 2 (Gray-Coded QPSK) ---
    rxSym_QPSK_Gray     = txSym_QPSK_Gray + noise_QPSK;
    det_b0_Gray         = real(rxSym_QPSK_Gray) > 0;
    det_b1_Gray         = imag(rxSym_QPSK_Gray) > 0;
    rxBits_QPSK_Gray    = reshape([det_b0_Gray, det_b1_Gray].', 1, []);
    
    BER_sim_QPSK_Gray(idx) = sum(rxBits_QPSK_Gray ~= txBits) / numBits;
    
    % --- Case B: Figure 3 (Non-Gray-Coded QPSK) ---
    rxSym_QPSK_NoGray   = txSym_QPSK_NoGray + noise_QPSK;
    det_I_NoGray        = real(rxSym_QPSK_NoGray) > 0;
    det_Q_NoGray        = imag(rxSym_QPSK_NoGray) > 0;
    
    det_b0_NoGray       = det_I_NoGray;
    det_b1_NoGray       = xor(det_I_NoGray, det_Q_NoGray);
    rxBits_QPSK_NoGray  = reshape([det_b0_NoGray, det_b1_NoGray].', 1, []);
    
    BER_sim_QPSK_NoGray(idx) = sum(rxBits_QPSK_NoGray ~= txBits) / numBits;
    
    % ---------------------------------------------------------------------
    % 3.3 8PSK Channel & Minimum-Distance Demapper (Eb = 1/3)
    % ---------------------------------------------------------------------
    Eb_8PSK             = 1 / 3;
    No_8PSK             = Eb_8PSK / current_EbNo;
    noise_8PSK          = sqrt(No_8PSK / 2) * ...
                          (randn(size(txSym_8PSK)) + 1j * randn(size(txSym_8PSK)));
    
    rxSym_8PSK          = txSym_8PSK + noise_8PSK;
    
    % Compute squared Euclidean distance to all 8 reference constellation points
    euclideanDistSq     = abs(rxSym_8PSK - constellation_8PSK.').^2;
    [~, minDistIdx]     = min(euclideanDistSq, [], 2);
    decVal_8PSK         = minDistIdx - 1;                % Convert to [0..7]
    
    % Recover binary triplet [b0, b1, b2] from decimal symbol index
    det_b0_8PSK         = floor(decVal_8PSK / 4);
    det_b1_8PSK         = floor(mod(decVal_8PSK, 4) / 2);
    det_b2_8PSK         = mod(decVal_8PSK, 2);
    rxBits_8PSK         = reshape([det_b0_8PSK, det_b1_8PSK, det_b2_8PSK].', 1, []);
    
    BER_sim_8PSK(idx)   = sum(rxBits_8PSK ~= txBits) / numBits;
    
    % ---------------------------------------------------------------------
    % 3.4 16-QAM Channel & Decision-Boundary Demapper (Eb = 2.5)
    % ---------------------------------------------------------------------
    Eb_16QAM            = 2.5;
    No_16QAM            = Eb_16QAM / current_EbNo;
    noise_16QAM         = sqrt(No_16QAM / 2) * ...
                          (randn(size(txSym_16QAM)) + 1j * randn(size(txSym_16QAM)));
    
    rxSym_16QAM         = txSym_16QAM + noise_16QAM;
    rx_I                = real(rxSym_16QAM);
    rx_Q                = imag(rxSym_16QAM);
    
    % Optimal ML decision boundaries for {-3, -1, +1, +3} are at 0 and +/-2
    det_b0_16QAM        = (rx_I >= 0);
    det_b1_16QAM        = (abs(rx_I) < 2);
    det_b2_16QAM        = (rx_Q >= 0);
    det_b3_16QAM        = (abs(rx_Q) < 2);
    
    rxBits_16QAM        = reshape([det_b0_16QAM, det_b1_16QAM, ...
                                   det_b2_16QAM, det_b3_16QAM].', 1, []);
                               
    BER_sim_16QAM(idx)  = sum(rxBits_16QAM ~= txBits) / numBits;
    
    % ---------------------------------------------------------------------
    % 3.5 BFSK Channel & Correlator Decision Demapper (Task 1.5, Eb = 1)
    % ---------------------------------------------------------------------
    No_BFSK             = Eb_BFSK / current_EbNo;
    noise_BFSK          = sqrt(No_BFSK / 2) * ...
                          (randn(size(txSym_BFSK)) + 1j * randn(size(txSym_BFSK)));
    
    rxSym_BFSK          = txSym_BFSK + noise_BFSK;
    
    % Optimal ML decision rule: compare projection on phi_2 vs. phi_1
    r1_BFSK             = real(rxSym_BFSK);
    r2_BFSK             = imag(rxSym_BFSK);
    rxBits_BFSK         = (r2_BFSK > r1_BFSK);
    
    BER_sim_BFSK(idx)   = sum(rxBits_BFSK ~= txBits) / numBits;
    
end


%% ========================================================================
%  4. THEORETICAL BER EXPRESSIONS
%  ========================================================================
%  Note: Using Q(x) = 0.5 * erfc(x / sqrt(2))

% 4.1 BPSK Exact Theoretical BER
BER_theo_BPSK           = 0.5 * erfc(sqrt(EbNo_lin));

% 4.2 QPSK (Gray Coded - Fig 2) Exact Theoretical BER
BER_theo_QPSK_Gray      = 0.5 * erfc(sqrt(EbNo_lin));

% 4.3 QPSK (Non-Gray Coded - Fig 3) Exact Theoretical BER
%     Horizontal symbol errors cause 2 bit errors; vertical cause 1 bit error
p_branch                = 0.5 * erfc(sqrt(EbNo_lin));
BER_theo_QPSK_NoGray    = 1.5 * p_branch - (p_branch .^ 2);

% 4.4 8PSK (Gray Coded) Tight Upper Bound Approximation
BER_theo_8PSK           = (1 / 3) * erfc(sqrt(3 * EbNo_lin) * sin(pi / 8));

% 4.5 16-QAM (Gray Coded) Tight Theoretical Approximation
BER_theo_16QAM          = (3 / 8) * erfc(sqrt((2 / 5) * EbNo_lin));

% 4.6 BFSK (Coherent Orthogonal) Exact Theoretical BER (Task 1.5)
BER_theo_BFSK           = 0.5 * erfc(sqrt(EbNo_lin / 2));


%% ========================================================================
%  5. TASK 1.5: BFSK BASEBAND EQUIVALENT WAVEFORM & PSD SIMULATION
%  ========================================================================
%  Passband Signal:
%    s_i(t) = sqrt(2*Eb/Tb) * cos(2*pi*f_i*t),  0 <= t <= Tb
%    f_i    = (n_c + i) / Tb,  i = 1, 2   =>  df = f_2 - f_1 = 1 / Tb
%
%  Baseband Equivalent Models:
%    Model A (Reference Carrier fc = f1 = (nc + 1)/Tb):
%      s_bb1(t) = sqrt(2*Eb/Tb)
%      s_bb2(t) = sqrt(2*Eb/Tb) * exp(j * 2*pi * t / Tb)
%
%    Model B (Symmetric Center Carrier fc = (f1 + f2)/2 = (nc + 1.5)/Tb):
%      s_bb1(t) = sqrt(2*Eb/Tb) * exp(-j * pi * t / Tb)
%      s_bb2(t) = sqrt(2*Eb/Tb) * exp(+j * pi * t / Tb)
%  ========================================================================

fprintf('Simulating Task 1.5 Baseband Equivalent BFSK Waveform & PSD...\n');

% Time & Sampling Parameters for Waveform Generation
Tb                      = 1.0;                      % Normalized bit period (1 s)
Fs                      = 64 / Tb;                  % Sampling frequency (64 samples/bit)
Ts                      = 1 / Fs;                   % Sampling interval
numBits_PSD             = 20000;                    % Number of bits for Welch PSD
t_bit                   = 0 : Ts : (Tb - Ts);       % Single-bit time vector [1 x 64]

% Extract bit subset and define complex envelope amplitude
bits_PSD                = txBits(1 : numBits_PSD).';
A_env                   = sqrt(2 * Eb_BFSK / Tb);

% --- 5.1 Model A Waveform: Carrier fc = f1 ---
s_bb1_ModelA            = A_env * ones(1, length(t_bit));
s_bb2_ModelA            = A_env * exp(1j * 2 * pi * (1 / Tb) * t_bit);

matWave_ModelA          = (bits_PSD == 0) * s_bb1_ModelA + ...
                          (bits_PSD == 1) * s_bb2_ModelA;
txWave_BFSK_ModelA      = reshape(matWave_ModelA.', 1, []);

% --- 5.2 Model B Waveform: Symmetric Carrier fc = (f1 + f2) / 2 ---
s_bb1_ModelB            = A_env * exp(-1j * pi * (1 / Tb) * t_bit);
s_bb2_ModelB            = A_env * exp(+1j * pi * (1 / Tb) * t_bit);

matWave_ModelB          = (bits_PSD == 0) * s_bb1_ModelB + ...
                          (bits_PSD == 1) * s_bb2_ModelB;
txWave_BFSK_ModelB      = reshape(matWave_ModelB.', 1, []);

% --- 5.3 Estimate Baseband PSD via Welch's Method ---
windowLen               = 2048;
nOverlap                = windowLen / 2;
nFFT                    = 4096;

[psd_ModelA, f_vec]     = pwelch(txWave_BFSK_ModelA, hamming(windowLen), ...
                                 nOverlap, nFFT, Fs, 'centered');
[psd_ModelB, ~]         = pwelch(txWave_BFSK_ModelB, hamming(windowLen), ...
                                 nOverlap, nFFT, Fs, 'centered');

% Normalized frequency vector (f * Tb)
f_norm                  = f_vec * Tb;


%% ========================================================================
%  6. DISPLAY SUMMARY TABLE IN COMMAND WINDOW
%  ========================================================================

fprintf('\n-------------------------------------------------------------------------------------\n');
fprintf(' Eb/No (dB) |    BPSK    | QPSK (Fig2)| QPSK (Fig3)|   16-QAM   |    8PSK    |    BFSK    \n');
fprintf('-------------------------------------------------------------------------------------\n');
for idx = 1 : numSNR
    fprintf('   %5.1f    |  %.3e |  %.3e |  %.3e |  %.3e |  %.3e |  %.3e \n', ...
        EbNo_dB(idx), BER_sim_BPSK(idx), BER_sim_QPSK_Gray(idx), ...
        BER_sim_QPSK_NoGray(idx), BER_sim_16QAM(idx), BER_sim_8PSK(idx), BER_sim_BFSK(idx));
end
fprintf('-------------------------------------------------------------------------------------\n\n');


%% ========================================================================
%  7. PLOTTING & VISUALIZATION (TASKS 1.4 & 1.5)
%  ========================================================================

% Common plot styling parameters
axisFontSize            = 11;
titleFontSize           = 12;
lineWidthTheo           = 1.8;
lineWidthSim            = 1.4;
markerSizeSim           = 7;
yAxisLimits             = [1e-6, 1];
xAxisLimits             = [EbNo_dB(1), EbNo_dB(end)];

% -------------------------------------------------------------------------
% Figure 1: BPSK Performance (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 1: BPSK BER Performance', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_BPSK, 'b-',  'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_BPSK,  'ro',  'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
grid on; box on;
xlim(xAxisLimits); ylim(yAxisLimits);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER vs. $E_b/N_0$ Performance for BPSK Modulation', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Theoretical BPSK', 'Simulated BPSK'}, 'Location', 'southwest', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 2: QPSK (Figure 2 - Gray Coded) Performance (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 2: QPSK (Gray) BER Performance', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_QPSK_Gray, 'r-', 'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_QPSK_Gray,  'ks', 'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
grid on; box on;
xlim(xAxisLimits); ylim(yAxisLimits);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER vs. $E_b/N_0$ Performance for Gray-Coded QPSK (Figure 2)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Theoretical QPSK (Gray)', 'Simulated QPSK (Gray)'}, 'Location', 'southwest', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 3: 8PSK (Gray Coded) Performance (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 3: 8PSK BER Performance', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_8PSK, 'm-', 'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_8PSK,  'kd', 'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
grid on; box on;
xlim(xAxisLimits); ylim(yAxisLimits);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER vs. $E_b/N_0$ Performance for Gray-Coded 8PSK', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Theoretical 8PSK (Tight Bound)', 'Simulated 8PSK'}, 'Location', 'southwest', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 4: 16-QAM (Gray Coded) Performance (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 4: 16-QAM BER Performance', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_16QAM, 'g-', 'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_16QAM,  'k^', 'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
grid on; box on;
xlim(xAxisLimits); ylim(yAxisLimits);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER vs. $E_b/N_0$ Performance for Gray-Coded 16-QAM', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Theoretical 16-QAM', 'Simulated 16-QAM'}, 'Location', 'southwest', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 5: Combined Comparison of All 4 Modulation Schemes (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 5: All Modulation Schemes Combined', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_BPSK,      'b-',  'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_BPSK,       'bo',  'LineWidth', lineWidthSim, 'MarkerSize', 6);
semilogy(EbNo_dB, BER_theo_QPSK_Gray, 'r--', 'LineWidth', lineWidthTheo);
semilogy(EbNo_dB, BER_sim_QPSK_Gray,  'rx',  'LineWidth', lineWidthSim, 'MarkerSize', 8);
semilogy(EbNo_dB, BER_theo_16QAM,     'g-',  'LineWidth', lineWidthTheo);
semilogy(EbNo_dB, BER_sim_16QAM,      'gs',  'LineWidth', lineWidthSim, 'MarkerSize', 6);
semilogy(EbNo_dB, BER_theo_8PSK,      'm-',  'LineWidth', lineWidthTheo);
semilogy(EbNo_dB, BER_sim_8PSK,       'md',  'LineWidth', lineWidthSim, 'MarkerSize', 6);
grid on; box on;
xlim(xAxisLimits); ylim([1e-5, 1]);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER Performance Comparison: BPSK, QPSK, 16-QAM, and 8PSK', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'BPSK (Theoretical)',   'BPSK (Simulated)', ...
        'QPSK (Theoretical)',   'QPSK (Simulated)', ...
        '16-QAM (Theoretical)', '16-QAM (Simulated)', ...
        '8PSK (Theoretical)',   '8PSK (Simulated)'}, ...
       'Location', 'southwest', 'FontSize', 9);

% -------------------------------------------------------------------------
% Figure 6: QPSK Comparison - Figure 2 (Gray) vs. Figure 3 (Non-Gray) (Task 1.4)
% -------------------------------------------------------------------------
figure('Name', 'Fig 6: QPSK Gray vs Non-Gray', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_QPSK_Gray,   'b-',  'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_QPSK_Gray,    'bo',  'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
semilogy(EbNo_dB, BER_theo_QPSK_NoGray, 'r--', 'LineWidth', lineWidthTheo);
semilogy(EbNo_dB, BER_sim_QPSK_NoGray,  'rs',  'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
grid on; box on;
xlim(xAxisLimits); ylim([1e-5, 1]);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('Effect of Bit Mapping on QPSK: Gray (Fig. 2) vs. Non-Gray (Fig. 3)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'QPSK Fig. 2 (Gray) - Theoretical',     'QPSK Fig. 2 (Gray) - Simulated', ...
        'QPSK Fig. 3 (Non-Gray) - Theoretical', 'QPSK Fig. 3 (Non-Gray) - Simulated'}, ...
       'Location', 'southwest', 'FontSize', 9.5);

% -------------------------------------------------------------------------
% Figure 7: BFSK BER vs. Eb/No Performance (Task 1.5)
% -------------------------------------------------------------------------
figure('Name', 'Fig 7: BFSK BER Performance (Task 1.5)', 'NumberTitle', 'off', 'Color', 'w');
semilogy(EbNo_dB, BER_theo_BFSK, 'b-',  'LineWidth', lineWidthTheo); hold on;
semilogy(EbNo_dB, BER_sim_BFSK,  'ro',  'LineWidth', lineWidthSim, 'MarkerSize', markerSizeSim);
semilogy(EbNo_dB, BER_theo_BPSK, 'k--', 'LineWidth', 1.4);
grid on; box on;
xlim(xAxisLimits); ylim(yAxisLimits);
set(gca, 'FontSize', axisFontSize);
xlabel('$E_b/N_0\ \mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Bit Error Rate (BER)', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('BER vs. $E_b/N_0$ Performance for Coherent Orthogonal BFSK', 'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Theoretical BFSK', 'Simulated BFSK', 'Theoretical BPSK (3 dB Reference)'}, ...
       'Location', 'southwest', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 8: BFSK Baseband Equivalent PSD - Carrier fc = f1 (Task 1.5)
% -------------------------------------------------------------------------
figure('Name', 'Fig 8: BFSK Baseband PSD (fc = f1)', 'NumberTitle', 'off', 'Color', 'w');
plot(f_norm, 10 * log10(psd_ModelA / max(psd_ModelA)), 'b-', 'LineWidth', 1.5);
grid on; box on;
xlim([-4, 4]); ylim([-50, 5]);
set(gca, 'FontSize', axisFontSize);
xlabel('Normalized Frequency $(f - f_1)T_b$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Normalized PSD $\mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('Simulated Baseband PSD of BFSK ($f_c = f_1$, Tones at $0$ and $1/T_b$)', ...
      'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Simulated Baseband PSD ($f_c = f_1$)'}, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);

% -------------------------------------------------------------------------
% Figure 9: BFSK Baseband Equivalent PSD - Symmetric Carrier fc = (f1+f2)/2 (Task 1.5)
% -------------------------------------------------------------------------
figure('Name', 'Fig 9: BFSK Baseband PSD (Symmetric fc)', 'NumberTitle', 'off', 'Color', 'w');
plot(f_norm, 10 * log10(psd_ModelB / max(psd_ModelB)), 'r-', 'LineWidth', 1.5);
grid on; box on;
xlim([-4, 4]); ylim([-50, 5]);
set(gca, 'FontSize', axisFontSize);
xlabel('Normalized Frequency $(f - f_c)T_b$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
ylabel('Normalized PSD $\mathrm{(dB)}$', 'Interpreter', 'latex', 'FontSize', titleFontSize);
title('Simulated Baseband PSD of BFSK ($f_c = \frac{f_1+f_2}{2}$, Tones at $\pm \frac{1}{2T_b}$)', ...
      'Interpreter', 'latex', 'FontSize', titleFontSize);
legend({'Simulated Baseband PSD ($f_c = (f_1+f_2)/2$)'}, 'Interpreter', 'latex', 'Location', 'northeast', 'FontSize', 10);