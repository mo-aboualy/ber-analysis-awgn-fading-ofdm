function results = ofdm_project_complete()
% OFDM_PROJECT_COMPLETE  End-to-end OFDM link simulation with diagnostics.
%
%   results = OFDM_PROJECT_COMPLETE() simulates an OFDM link and returns a
%   structure containing BER, EVM, error counts, and diagnostic snapshots.
%
%   Configuration
%     Modulation : BPSK and Gray-coded 16-QAM
%     Channels   : (1) flat Rayleigh fading (one complex gain per OFDM symbol)
%                  (2) frequency-selective Rayleigh fading, modeled as an
%                      independent complex Rayleigh gain on every subcarrier
%                      (frequency-domain model specified in the assignment)
%     Coding     : uncoded, and rate-1/5 repetition code (majority decoding)
%     Receiver   : hard-decision demapping with perfect channel knowledge (CSI)
%
%   Outputs
%     - BER vs. SNR figures (one per modulation, four curves each)
%     - EVM vs. SNR figures (one per modulation, four curves each)
%     - Transmitter/receiver diagnostic figures per modulation-channel pair
%     - BER heatmap covering all 8 scenarios
%     - Results saved to 'ofdm_project_results.mat'; every figure saved as PNG
%
%   Implementation notes
%     - Only base MATLAB functions (fft/ifft, plotting) are used; the
%       mapper/demapper are implemented manually (no qammod/qamdemod/ofdmmod).
%     - The OFDM grid is a time-frequency map: rows are subcarriers
%       (frequency), columns are OFDM symbols (time), and each cell holds
%       one modulation symbol.
%     - Transmit chain: bits -> coding -> interleaving -> mapping -> IFFT ->
%       cyclic prefix -> channel. The receiver applies the inverse steps.
%     - The cyclic prefix (CP) copies the tail of each useful symbol to its
%       front. When the CP is at least as long as the channel's excess delay,
%       it prevents inter-symbol interference and preserves subcarrier
%       orthogonality.
%     - Pilot-based estimation (H = Ypilot/Xpilot, interpolated over the grid)
%       is not used for the BER curves. The assignment specifies independently
%       faded subcarriers, which sparse-pilot interpolation cannot track, so
%       perfect CSI is used as the baseline.
%     - The assignment requires a 256-point IFFT with BER curves for BPSK and
%       16-QAM (4 cases each: 2 channels x coded/uncoded).

clc; close all; rng(42, 'twister');

% Output folder: every figure is saved as an individual PNG file.
figDir = fullfile(pwd, 'OFDM_figures');
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

%% 1. Simulation parameters
Nfft   = 256;                  % IFFT/FFT size (assignment requirement)
cpLen  = 16;                   % CP length in samples (assumption: not specified)
nOFDMSymbolsPerFrame = 8;      % OFDM symbols per simulated frame
nFrames = 100;                 % Frames per SNR point (use 500-1000 for final curves)
snrDb  = 0:2:24;               % Simulated SNR range [dB]
diagnosticSNRdB = 12;          % SNR at which diagnostic snapshots are captured

MList        = [2 16];
modNames     = {'BPSK','16-QAM'};
channelNames = {'Flat Rayleigh fading', ...
                'Frequency-selective Rayleigh fading'};
codingNames  = {'No coding','Rate-1/5 repetition'};

% Block interleaver dimensions [rows cols]; rows*cols must equal the number
% of coded bits per OFDM symbol, Nfft*log2(M). 16-QAM (32x32) follows the
% assignment; BPSK is unspecified, so 16x16 (256 bits) is assumed.
interleaverSizes = [16 16; 32 32];

assert(cpLen >= 0 && cpLen < Nfft, 'cpLen must be between 0 and Nfft-1.');
assert(all(interleaverSizes(:,1).*interleaverSizes(:,2) == ...
    Nfft*log2(MList(:))), ...
    'Interleaver size must equal the number of bits per OFDM symbol.');

nMod  = numel(MList);
nChan = numel(channelNames);
nCode = numel(codingNames);
nSNR  = numel(snrDb);

% Result arrays are indexed as (modulation, channel, coding, SNR).
BER         = zeros(nMod,nChan,nCode,nSNR);
EVM         = zeros(nMod,nChan,nCode,nSNR);
errorCounts = zeros(size(BER));
bitCounts   = zeros(size(BER));

diagSNRidx = find(snrDb == diagnosticSNRdB,1);
if isempty(diagSNRidx)
    [~,diagSNRidx] = min(abs(snrDb-diagnosticSNRdB));
end

% One uncoded diagnostic snapshot per modulation/channel pair.
diagnostics = cell(nMod,nChan);

%% 2. Monte Carlo simulation
for im = 1:nMod
    M = MList(im);
    interRows = interleaverSizes(im,1);
    interCols = interleaverSizes(im,2);

    fprintf('\n=== %s ===\n', modNames{im});

    for ich = 1:nChan
        for ic = 1:nCode
            useCoding = (ic == 2);
            fprintf('  %-36s | %-22s\n', channelNames{ich}, codingNames{ic});

            for isnr = 1:nSNR
                totalErrors = 0;
                totalBits   = 0;
                evmSum      = 0;

                for frameIdx = 1:nFrames
                    % Keep a snapshot from the first uncoded frame at the
                    % diagnostic SNR.
                    keepDiag = (ic == 1 && isnr == diagSNRidx && frameIdx == 1);

                    [nErr,nBits,evmFrame,diagFrame] = simulateOneFrame( ...
                        M,useCoding,ich,Nfft,cpLen,nOFDMSymbolsPerFrame, ...
                        10^(snrDb(isnr)/10),interRows,interCols,keepDiag);

                    totalErrors = totalErrors + nErr;
                    totalBits   = totalBits + nBits;
                    evmSum      = evmSum + evmFrame;

                    if keepDiag
                        diagnostics{im,ich} = diagFrame;
                    end
                end

                BER(im,ich,ic,isnr)         = totalErrors/max(totalBits,1);
                EVM(im,ich,ic,isnr)         = evmSum/nFrames;
                errorCounts(im,ich,ic,isnr) = totalErrors;
                bitCounts(im,ich,ic,isnr)   = totalBits;
            end
        end
    end
end

%% 3. Diagnostic figures (transmitter, channel, receiver)
% Each figure holds at most three axes and is saved as its own PNG.
for im = 1:nMod
    for ich = 1:nChan
        d = diagnostics{im,ich};
        if ~isempty(d)
            plotDiagnostics(d, MList(im), modNames{im}, ...
                channelNames{ich}, diagnosticSNRdB, Nfft, cpLen, figDir);
        end
    end
end

%% 4. BER versus SNR (one figure per modulation, four curves each)
curveStyles = {'-o','-s','-^','-d'};
for im = 1:nMod
    createFigure(['BER vs SNR - ' modNames{im}], 900, 600);
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    legendText = cell(1,4);
    curveIdx = 0;

    for ich = 1:nChan
        for ic = 1:nCode
            curveIdx = curveIdx + 1;
            y = reshape(BER(im,ich,ic,:),1,[]);
            y(y == 0) = NaN;   % Zero-error points cannot be shown on a log axis.
            semilogy(ax,snrDb,y,curveStyles{curveIdx}, ...
                'LineWidth',1.6,'MarkerSize',6);
            legendText{curveIdx} = [channelNames{ich} ', ' codingNames{ic}];
        end
    end

    set(ax,'YScale','log','FontSize',11);
    xlabel(ax,'SNR (dB)');
    ylabel(ax,'Bit Error Rate (BER)');
    title(ax,['BER versus SNR - ' modNames{im}]);
    legend(ax,legendText,'Location','southwest','FontSize',9);
    xlim(ax,[min(snrDb) max(snrDb)]);
    ylim(ax,[1e-5 1]);
    hold(ax,'off');
    saveCurrentFigure(figDir, ['BER_vs_SNR_' slugify(modNames{im})]);
end

%% 5. EVM versus SNR (one figure per modulation, four curves each)
for im = 1:nMod
    createFigure(['EVM vs SNR - ' modNames{im}], 900, 600);
    ax = axes; hold(ax,'on'); grid(ax,'on'); box(ax,'on');
    legendText = cell(1,4);
    curveIdx = 0;

    for ich = 1:nChan
        for ic = 1:nCode
            curveIdx = curveIdx + 1;
            y = reshape(EVM(im,ich,ic,:),1,[]);
            plot(ax,snrDb,y,curveStyles{curveIdx}, ...
                'LineWidth',1.6,'MarkerSize',6);
            legendText{curveIdx} = [channelNames{ich} ', ' codingNames{ic}];
        end
    end

    set(ax,'FontSize',11);
    xlabel(ax,'SNR (dB)');
    ylabel(ax,'RMS EVM (%)');
    title(ax,['RMS EVM versus SNR - ' modNames{im}]);
    legend(ax,legendText,'Location','northeast','FontSize',9);
    xlim(ax,[min(snrDb) max(snrDb)]);
    hold(ax,'off');
    saveCurrentFigure(figDir, ['EVM_vs_SNR_' slugify(modNames{im})]);
end

%% 6. BER heatmap across all scenarios
nScen = nMod*nChan*nCode;
scenarioNames = cell(nScen,1);
berMatrix     = zeros(nScen,nSNR);
row = 0;
for im = 1:nMod
    for ich = 1:nChan
        for ic = 1:nCode
            row = row + 1;
            scenarioNames{row} = [modNames{im} ' | ' ...
                channelNames{ich} ' | ' codingNames{ic}];
            berMatrix(row,:) = reshape(BER(im,ich,ic,:),1,[]);
        end
    end
end

createFigure('BER Heatmap - All Scenarios', 1200, 600);
imagesc(snrDb,1:size(berMatrix,1),log10(max(berMatrix,1e-6)));
axis xy; colorbar;
xlabel('SNR (dB)');
ylabel('Scenario');
title(sprintf('log_{10}(BER) Heatmap - All %d Scenarios',nScen));
yticks(1:numel(scenarioNames));
yticklabels(scenarioNames);
set(gca,'TickLabelInterpreter','none','FontSize',8);
saveCurrentFigure(figDir, sprintf('BER_heatmap_all_%d_scenarios',nScen));

%% 7. Console summary and results export
fprintf('\n================ BER SUMMARY (SNR = %g dB) ================\n', ...
    snrDb(end));
for im = 1:nMod
    fprintf('\n%s\n',modNames{im});
    for ich = 1:nChan
        for ic = 1:nCode
            curve = reshape(BER(im,ich,ic,:),1,[]);
            fprintf('  %-36s | %-22s | BER = %.6g\n', ...
                channelNames{ich},codingNames{ic},curve(end));
        end
    end
end

results.BER = BER;
results.EVM = EVM;
results.errorCounts = errorCounts;
results.bitCounts = bitCounts;
results.snrDb = snrDb;
results.MList = MList;
results.modNames = modNames;
results.channelNames = channelNames;
results.codingNames = codingNames;
results.Nfft = Nfft;
results.cpLen = cpLen;
results.nFrames = nFrames;
results.nOFDMSymbolsPerFrame = nOFDMSymbolsPerFrame;
results.diagnostics = diagnostics;
save('ofdm_project_results.mat','results');

fprintf('\nResults saved to ofdm_project_results.mat\n');
fprintf('Figures saved in: %s\n', figDir);
end

%% ========================================================================
function [nErrors,nInfoBits,evmPercent,d] = simulateOneFrame( ...
    M,useCoding,channelType,Nfft,cpLen,nOFDMSymbols,snrLinear, ...
    interRows,interCols,keepDiagnostics)
% SIMULATEONEFRAME  Simulate one OFDM frame with perfect-CSI equalization.
%
%   channelType: 1 = flat Rayleigh, 2 = frequency-selective Rayleigh.
%   Returns the bit-error count, number of information bits, RMS EVM [%],
%   and (if keepDiagnostics is true) a structure of intermediate signals.

k = log2(M);
bitsPerOFDMSymbol = Nfft*k;
if useCoding
    nInfoBitsPerSymbol = floor(bitsPerOFDMSymbol/5);   % rate-1/5 repetition
else
    nInfoBitsPerSymbol = bitsPerOFDMSymbol;
end

infoBits = false(nInfoBitsPerSymbol,nOFDMSymbols);
txGrid   = complex(zeros(Nfft,nOFDMSymbols));

% ---- Transmitter: bits -> coding -> interleaving -> symbol mapping ----
for s = 1:nOFDMSymbols
    b = logical(randi([0 1],nInfoBitsPerSymbol,1));
    infoBits(:,s) = b;

    if useCoding
        % Repeat each information bit 5 times and zero-pad to the full
        % symbol capacity (e.g., BPSK: 51 bits -> 255 coded + 1 zero = 256).
        encoded = repelem(double(b),5);
        encoded = [encoded; zeros(bitsPerOFDMSymbol-numel(encoded),1)];
    else
        encoded = double(b);
    end

    encodedInterleaved = blockInterleave(encoded,interRows,interCols);
    txGrid(:,s) = mapBitsToSymbols(encodedInterleaved,M);
end

% ---- IFFT and cyclic prefix ----
% ifftshift converts the centered subcarrier ordering to MATLAB's FFT
% ordering. The sqrt(Nfft) factor keeps the time-domain power at unity.
txUseful = ifft(ifftshift(txGrid,1),Nfft,1)*sqrt(Nfft);
if cpLen > 0
    % Cyclic prefix: copy the last cpLen samples to the front of the symbol.
    txBlocks = [txUseful(end-cpLen+1:end,:); txUseful];
else
    txBlocks = txUseful;
end

% Noise variance is set from the measured transmit power and the target SNR.
signalPower = mean(abs(txBlocks(:)).^2);
noiseVar    = signalPower/snrLinear;

rxBlocks = complex(zeros(size(txBlocks)));
rxGrid   = complex(zeros(size(txGrid)));
Htrue    = complex(zeros(size(txGrid)));

% ---- Channel, noise, CP removal, and FFT (per OFDM symbol) ----
for s = 1:nOFDMSymbols
    if channelType == 1
        % Flat Rayleigh fading: one complex gain shared by all subcarriers
        % of this OFDM symbol; independent from symbol to symbol.
        h = (randn+1i*randn)/sqrt(2);
        Htrue(:,s) = h;
        channelOut = h*txBlocks(:,s);
    else
        % Frequency-selective Rayleigh fading (assignment model): each
        % subcarrier k has an independent complex Rayleigh gain H[k],
        % applied directly in the frequency domain.
        H = (randn(Nfft,1)+1i*randn(Nfft,1))/sqrt(2);
        Htrue(:,s) = H;
        Yfreq   = H.*txGrid(:,s);
        yUseful = ifft(ifftshift(Yfreq,1),Nfft,1)*sqrt(Nfft);
        if cpLen > 0
            channelOut = [yUseful(end-cpLen+1:end); yUseful];
        else
            channelOut = yUseful;
        end
    end

    % Additive white complex Gaussian noise.
    noise = sqrt(noiseVar/2)*(randn(size(channelOut))+1i*randn(size(channelOut)));
    rxBlocks(:,s) = channelOut + noise;

    % Discard the CP, then transform back to the subcarrier domain.
    if cpLen > 0
        yUsefulRx = rxBlocks(cpLen+1:cpLen+Nfft,s);
    else
        yUsefulRx = rxBlocks(:,s);
    end
    rxGrid(:,s) = fftshift(fft(yUsefulRx,Nfft,1),1)/sqrt(Nfft);
end

% ---- Equalization (perfect CSI) ----
% Zero-forcing equalizer written as Y*conj(H)/(|H|^2 + eps). The small
% regularization term only prevents division by zero for deep fades.
regularization = 1e-10;
eqGrid = rxGrid.*conj(Htrue)./(abs(Htrue).^2+regularization);

% RMS EVM of the equalized symbols relative to the transmitted symbols.
evmPercent = sqrt(mean(abs(eqGrid(:)-txGrid(:)).^2) / ...
    max(mean(abs(txGrid(:)).^2),eps))*100;

% ---- Receiver: demapping -> deinterleaving -> (optional) decoding ----
% Repetition decoding uses a majority vote over each group of 5 bits.
nErrors   = 0;
nInfoBits = 0;
for s = 1:nOFDMSymbols
    codedRx = demapSymbolsToBits(eqGrid(:,s),M);
    decodedInterleaver = blockDeinterleave(codedRx,interRows,interCols);

    if useCoding
        requiredCodedBits = 5*nInfoBitsPerSymbol;
        groups = reshape(decodedInterleaver(1:requiredCodedBits),5,[]);
        recoveredBits = sum(groups,1).' >= 3;
    else
        recoveredBits = decodedInterleaver;
    end

    original      = infoBits(:,s);
    recoveredBits = logical(recoveredBits(:));
    nErrors   = nErrors + sum(recoveredBits ~= original);
    nInfoBits = nInfoBits + numel(original);
end

% ---- Optional snapshot of intermediate signals for diagnostic plots ----
d = [];
if keepDiagnostics
    d.txGrid     = txGrid;
    d.txBlocks   = txBlocks;
    d.rxBlocks   = rxBlocks;
    d.rxGrid     = rxGrid;
    d.eqGrid     = eqGrid;
    d.Htrue      = Htrue;
    d.evmPercent = evmPercent;
end
end

%% ========================================================================
function plotDiagnostics(d,M,modName,channelName,snrDb,Nfft,cpLen,figDir)
% PLOTDIAGNOSTICS  Generate four diagnostic figures for one scenario.
%
%   A: transmitter       - ideal constellation, OFDM grid, TX waveform
%   B: transmit/channel  - transmit spectrum, received waveform
%   C: channel response  - magnitude and phase in ONE image (2 plots)
%   D: receiver          - constellation before/after equalization
%
%   No figure contains more than three axes; each is saved as its own PNG.

idealLabels = (0:M-1).';
idealBits   = labelsToBits(idealLabels,log2(M));
idealConstellation = mapBitsToSymbols(idealBits,M);
tag = [slugify(modName) '_' slugify(channelName)];
nSamples = min(size(d.txBlocks,1),Nfft+cpLen);

%% Figure A: ideal constellation, transmitted grid, TX waveform
createFigure([modName ' diagnostics A - Transmitter'], 1400, 450);

subplot(1,3,1);
plot(real(idealConstellation),imag(idealConstellation),'bo', ...
    'MarkerFaceColor','b','MarkerSize',6);
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title(['Ideal constellation: ' modName]);

subplot(1,3,2);
imagesc(1:size(d.txGrid,2),1:Nfft,abs(d.txGrid));
axis xy; colorbar;
xlabel('OFDM symbol index'); ylabel('Subcarrier index');
title('Transmitted OFDM grid |X[k,m]|');

subplot(1,3,3);
plot(0:nSamples-1,real(d.txBlocks(1:nSamples,1)),'b-'); hold on;
plot(0:nSamples-1,imag(d.txBlocks(1:nSamples,1)),'r--'); hold off;
grid on; xlabel('Sample index'); ylabel('Amplitude');
title('TX waveform (first OFDM symbol)');
legend('Real','Imaginary','Location','best');

sgtitle([modName ' | ' channelName ' | Transmitter diagnostics | SNR = ' ...
    num2str(snrDb) ' dB']);
saveCurrentFigure(figDir,['Diagnostics_A_Transmitter_' tag]);

%% Figure B: transmit spectrum and received waveform
createFigure([modName ' diagnostics B - Transmit/Channel'], 1100, 450);

tx    = d.txBlocks(:);
nSpec = 4*2^nextpow2(numel(tx));
spec  = fftshift(fft(tx,nSpec));
sampleRate = Nfft*15000;       % Assumes 15 kHz subcarrier spacing (LTE-like)
freq   = (-nSpec/2:nSpec/2-1)*(sampleRate/nSpec);
specDb = 20*log10(abs(spec)/(max(abs(spec))+eps)+1e-12);

subplot(1,2,1);
plot(freq/1e6,specDb,'LineWidth',1);
grid on; xlabel('Frequency (MHz)'); ylabel('Normalized magnitude (dB)');
title('Transmit spectrum');

subplot(1,2,2);
plot(0:nSamples-1,real(d.rxBlocks(1:nSamples,1)),'b-'); hold on;
plot(0:nSamples-1,imag(d.rxBlocks(1:nSamples,1)),'r--'); hold off;
grid on; xlabel('Sample index'); ylabel('Amplitude');
title('RX waveform (first OFDM symbol)');
legend('Real','Imaginary','Location','best');

sgtitle([modName ' | ' channelName ' | Transmit/channel diagnostics | SNR = ' ...
    num2str(snrDb) ' dB']);
saveCurrentFigure(figDir,['Diagnostics_B_Transmit_Channel_' tag]);

%% Figure C: channel magnitude and phase in one image (two plots)
% The wrapped phase is plotted: with independent per-subcarrier fading the
% phase is uncorrelated across k, so unwrapping would be meaningless.
createFigure([modName ' channel magnitude and phase'], 1200, 450);
H = d.Htrue(:,1);

subplot(1,2,1);
plot(1:Nfft,abs(H),'b-','LineWidth',1.1);
grid on; xlim([1 Nfft]);
xlabel('Subcarrier index'); ylabel('|H[k]|');
title('Channel magnitude response');

subplot(1,2,2);
plot(1:Nfft,angle(H),'r.-','LineWidth',0.8,'MarkerSize',6);
grid on; xlim([1 Nfft]); ylim([-pi pi]);
yticks([-pi -pi/2 0 pi/2 pi]);
yticklabels({'-\pi','-\pi/2','0','\pi/2','\pi'});
xlabel('Subcarrier index'); ylabel('Phase (rad)');
title('Channel phase response');

sgtitle([modName ' | ' channelName ' | Channel magnitude and phase']);
saveCurrentFigure(figDir,['Channel_Magnitude_and_Phase_' tag]);

%% Figure D: constellations before and after equalization
createFigure([modName ' diagnostics D - Receiver'], 1100, 500);

subplot(1,2,1);
plot(real(d.rxGrid(:)),imag(d.rxGrid(:)),'.','Color',[0.85 0.25 0.2]);
axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title('Before equalization');

subplot(1,2,2);
plot(real(d.eqGrid(:)),imag(d.eqGrid(:)),'.','Color',[0.1 0.45 0.85]);
hold on;
plot(real(idealConstellation),imag(idealConstellation), ...
    'kx','MarkerSize',8,'LineWidth',1.4);
hold off; axis equal; grid on;
xlabel('In-phase'); ylabel('Quadrature');
title(['After equalization (EVM = ' num2str(d.evmPercent,'%.1f') '%)']);
legend('Equalized symbols','Ideal points','Location','best');

sgtitle([modName ' | ' channelName ' | Receiver constellation diagnostics']);
saveCurrentFigure(figDir,['Diagnostics_D_Receiver_' tag]);
end

%% ========================================================================
%  Utility functions
%  ========================================================================
function createFigure(figName,widthPx,heightPx)
% Create a white figure of a fixed pixel size so saved PNGs are consistent.
fig = figure('Name',figName,'Color','w','NumberTitle','off');
pos = fig.Position;
fig.Position = [pos(1) pos(2) widthPx heightPx];
end

function saveCurrentFigure(figDir,fileStem)
% Save the current figure as a PNG file in figDir.
if ~exist(figDir,'dir')
    mkdir(figDir);
end
drawnow;
saveas(gcf,fullfile(figDir,[fileStem '.png']));
end

function out = slugify(in)
% Convert a title into a lowercase, file-name-safe string.
out = lower(regexprep(char(in),'[^A-Za-z0-9]+','_'));
out = regexprep(out,'(^_+|_+$)','');
end

function bits = labelsToBits(labels,k)
% Convert integer labels to a bit stream (most significant bit first).
labels = labels(:);
bits = zeros(numel(labels)*k,1);
for r = 1:numel(labels)
    for b = 1:k
        bits((r-1)*k+b) = bitget(uint32(labels(r)),k-b+1);
    end
end
end

function symbols = mapBitsToSymbols(bits,M)
% Map bits to unit-average-power symbols: BPSK or Gray-coded 16-QAM.
bits = double(bits(:));
k = log2(M);
if mod(numel(bits),k) ~= 0
    error('Input bit count must be divisible by log2(M).');
end
B = reshape(bits,k,[]).';

switch M
    case 2
        symbols = 1 - 2*B(:,1);              % 0 -> +1, 1 -> -1
    case 16
        I = grayPairToLevel(B(:,1),B(:,2));
        Q = grayPairToLevel(B(:,3),B(:,4));
        symbols = (I + 1i*Q)/sqrt(10);
    otherwise
        error('Supported M values are 2 and 16.');
end
symbols = symbols(:);
end

function bits = demapSymbolsToBits(symbols,M)
% Hard-decision demapper; exact inverse of mapBitsToSymbols.
symbols = symbols(:);
switch M
    case 2
        bits = double(real(symbols) < 0);
    case 16
        BI = grayLevelToPair(real(symbols)*sqrt(10));
        BQ = grayLevelToPair(imag(symbols)*sqrt(10));
        bits = reshape([BI BQ].',[],1);
    otherwise
        error('Supported M values are 2 and 16.');
end
bits = logical(bits(:));
end

function level = grayPairToLevel(b1,b2)
% Gray mapping of a bit pair to a 4-PAM level: 00 -> -3, 01 -> -1, 11 -> +1, 10 -> +3.
level = zeros(size(b1));
level(b1==0 & b2==0) = -3;
level(b1==0 & b2==1) = -1;
level(b1==1 & b2==1) =  1;
level(b1==1 & b2==0) =  3;
end

function B = grayLevelToPair(x)
% Hard-decision inverse of grayPairToLevel (decision thresholds at -2, 0, +2).
B = zeros(numel(x),2);
idx = x < -2;          B(idx,:) = repmat([0 0],sum(idx),1);
idx = x >= -2 & x < 0; B(idx,:) = repmat([0 1],sum(idx),1);
idx = x >= 0 & x < 2;  B(idx,:) = repmat([1 1],sum(idx),1);
idx = x >= 2;          B(idx,:) = repmat([1 0],sum(idx),1);
end

function y = blockInterleave(x,nRows,nCols)
% Block interleaver: write column-wise, read row-wise. Spreads adjacent
% bits apart so burst errors are distributed before decoding.
A = reshape(x(:),nRows,nCols);
y = reshape(A.',[],1);
end

function x = blockDeinterleave(y,nRows,nCols)
% Inverse permutation of blockInterleave.
A = reshape(y(:),nCols,nRows).';
x = A(:);
end