# BER Analysis of Digital Modulation Schemes over AWGN, Rayleigh Flat-Fading, and Frequency-Selective Channels with OFDM

**A from-scratch MATLAB simulation study of digital modulation, from single-carrier AWGN links to Rayleigh fading and a full OFDM transceiver.**

Digital Communications · 3rd Year · Electronics & Communications Dept. · Faculty of Engineering, Cairo University

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-orange)
![Toolbox](https://img.shields.io/badge/requires-Signal%20Processing%20Toolbox-blue)
![Status](https://img.shields.io/badge/status-Part%201%20done%20%7C%20Parts%202--3%20in%20progress-yellow)

---

## Overview

This repository is a three-part simulation project that builds a digital communication link step by step. Every block (mapper, channel, demapper, BER counter) is written by hand in plain MATLAB with no `comm.*` toolbox objects, and every simulated curve is checked against a closed-form theoretical BER.

| Part | Topic | Channel | Status |
|------|-------|---------|--------|
| **1** | Single-carrier system: BPSK, QPSK (Gray and non-Gray), 8PSK, 16-QAM, BFSK, plus BFSK PSD | AWGN | Done |
| **2** | BER of BPSK, QPSK and 16-QAM | Rayleigh flat fading, channel known at the receiver | In progress |
| **3** | OFDM transceiver: coding, interleaver, mapper, 256-point IFFT, cyclic extension | Rayleigh flat fading and frequency-selective fading | In progress |

> **Deadline:** all parts due 14 Dec 2026, 11:59 PM.

---

## Table of Contents

1. [Part 1: Single-Carrier System](#part-1-single-carrier-system)
2. [Part 2: Rayleigh Flat Fading](#part-2-rayleigh-flat-fading) *(to be continued)*
3. [Part 3: OFDM System](#part-3-ofdm-system) *(to be continued)*
4. [Getting Started](#getting-started)
5. [Repository Structure](#repository-structure)
6. [Team](#team)

---

# Part 1: Single-Carrier System

Part 1 simulates the **baseband-equivalent** single-carrier chain (no carriers) over an **AWGN channel** and measures the bit-error rate (BER) versus $E_b/N_0$.

```mermaid
flowchart LR
    A[Random bits] --> B[Mapper]
    B --> C[AWGN channel]
    C --> D[Demapper]
    D --> E[BER calculator]
    A -.-> E
```

## What is implemented

| Scheme | Bits/symbol | Mapping | $E_s$ | $E_b$ | Demapper |
|--------|:-----------:|---------|:-----:|:-----:|----------|
| BPSK | 1 | $0 \to -1,\ 1 \to +1$ | 1 | 1 | sign of received sample |
| QPSK, Gray (Fig. 2) | 2 | $b_0 \to I,\ b_1 \to Q$ | 2 | 1 | independent sign decisions on $I$ and $Q$ |
| QPSK, non-Gray (Fig. 3) | 2 | $I = b_0,\ Q = b_0 \oplus b_1$ | 2 | 1 | sign decisions, then $b_1 = \hat I \oplus \hat Q$ |
| 8PSK | 3 | unit circle, phases from the Fig. 2 constellation | 1 | 1/3 | minimum Euclidean distance (ML) |
| 16-QAM, Gray | 4 | $\{00,01,11,10\} \to \{-3,-1,+1,+3\}$ per axis | 10 | 2.5 | per-axis thresholds at $0$ and $\pm 2$ (ML) |
| BFSK | 1 | orthonormal basis, $0 \to [\sqrt{E_b},0],\ 1 \to [0,\sqrt{E_b}]$ | 1 | 1 | correlator comparison $r_2 > r_1$ |

### Noise scaling

For every scheme the noise power is set from the target $E_b/N_0$ using that scheme's own $E_b$:

$$
N_0 = \frac{E_b}{E_b/N_0}, \qquad
n \sim \mathcal{N}\!\left(0, \tfrac{N_0}{2}\right) \text{ per real dimension}
$$

This keeps the comparison across modulation orders fair, since $E_b$, not $E_s$, is the common reference.

## Theoretical references

Simulated points are overlaid on these expressions, with $Q(x) = \tfrac{1}{2}\,\mathrm{erfc}\!\left(x/\sqrt{2}\right)$:

| Scheme | Theoretical BER | Type |
|--------|-----------------|------|
| BPSK | $\tfrac{1}{2}\,\mathrm{erfc}\!\left(\sqrt{E_b/N_0}\right)$ | exact |
| QPSK, Gray | $\tfrac{1}{2}\,\mathrm{erfc}\!\left(\sqrt{E_b/N_0}\right)$ | exact |
| QPSK, non-Gray | $\tfrac{3}{2}p - p^2,\quad p = \tfrac{1}{2}\,\mathrm{erfc}\!\left(\sqrt{E_b/N_0}\right)$ | exact |
| 8PSK, Gray | $\tfrac{1}{3}\,\mathrm{erfc}\!\left(\sqrt{3E_b/N_0}\,\sin\tfrac{\pi}{8}\right)$ | tight approximation |
| 16-QAM, Gray | $\tfrac{3}{8}\,\mathrm{erfc}\!\left(\sqrt{\tfrac{2}{5}\,E_b/N_0}\right)$ | tight approximation |
| BFSK, coherent orthogonal | $\tfrac{1}{2}\,\mathrm{erfc}\!\left(\sqrt{\tfrac{E_b}{2N_0}}\right)$ | exact |

**Why non-Gray QPSK is worse.** In Fig. 3, a decision error on the $I$ axis flips both bits (because $b_1$ depends on $I$), while an error on $Q$ alone flips one bit, and errors on both axes flip one bit net. Averaging over the two bits gives $\tfrac{3}{2}p - p^2$, which is about $1.5\times$ the Gray-coded BER at high SNR. Gray mapping ensures nearest-neighbour symbol errors cost exactly one bit error.

## BFSK (Task 1.5)

The transmitted tones are $f_i = (n_c + i)/T_b$ for $i = 1, 2$, so the tone spacing is $1/T_b$.

- **Basis functions:** $\phi_i(t) = \sqrt{2/T_b}\,\cos(2\pi f_i t)$, $0 \le t \le T_b$, which are orthonormal.
- **Signal space:** $s_1 = [\sqrt{E_b},\ 0]$ and $s_2 = [0,\ \sqrt{E_b}]$.
- **Baseband equivalents**, with $A = \sqrt{2E_b/T_b}$, for two choices of reference carrier:

| Model | Carrier $f_c$ | $s_{bb,1}(t)$ | $s_{bb,2}(t)$ | Tones at |
|-------|---------------|---------------|---------------|----------|
| A | $f_1$ | $A$ | $A\,e^{j2\pi t/T_b}$ | $0$ and $1/T_b$ |
| B | $(f_1+f_2)/2$ | $A\,e^{-j\pi t/T_b}$ | $A\,e^{+j\pi t/T_b}$ | $\pm 1/(2T_b)$ |

- **PSD:** estimated with Welch's method (`pwelch`, Hamming window of 2048 samples, 50% overlap, 4096-point FFT) on a 20,000-bit waveform sampled at 64 samples per bit, and plotted against normalized frequency $(f - f_c)T_b$ for both models.
- **Sanity check:** coherent orthogonal BFSK needs 3 dB more $E_b/N_0$ than BPSK for the same BER, so the BFSK figure includes the BPSK curve as a reference.

## Simulation parameters

| Parameter | Value | Note |
|-----------|-------|------|
| Bits per $E_b/N_0$ point | $1.2 \times 10^6$ | a multiple of 12, so it divides evenly by 1, 2, 3 and 4 bits/symbol with no padding |
| $E_b/N_0$ range | $-2$ to $14$ dB, 1 dB steps | |
| Random seed | `rng('default')` | fully reproducible |
| Same bit stream for all schemes | yes | differences between curves come from the scheme, not from the data |

## Output figures

Running the script produces nine figures:

| # | Figure | Task |
|---|--------|------|
| 1 | BPSK: simulated vs. theoretical | 1.4 |
| 2 | QPSK Gray (Fig. 2) | 1.4 |
| 3 | 8PSK Gray | 1.4 |
| 4 | 16-QAM Gray | 1.4 |
| 5 | BPSK, QPSK, 8PSK and 16-QAM on one graph | 1.4 |
| 6 | QPSK Gray (Fig. 2) vs. non-Gray (Fig. 3) | 1.4 |
| 7 | BFSK BER with BPSK reference | 1.5 |
| 8 | BFSK baseband PSD, $f_c = f_1$ | 1.5 |
| 9 | BFSK baseband PSD, $f_c = (f_1+f_2)/2$ | 1.5 |

A formatted table of all simulated BER values is also printed to the Command Window.

<!-- After running the script, export the figures into results/part1/ and uncomment:

![All schemes](results/part1/fig5_all_schemes.png)
![QPSK Gray vs non-Gray](results/part1/fig6_qpsk_gray_vs_nongray.png)
![BFSK PSD](results/part1/fig8_bfsk_psd.png)
-->

## Key observations

- **BPSK and Gray QPSK overlap exactly** on a BER vs. $E_b/N_0$ plot: QPSK is two independent BPSK streams on orthogonal axes with the same energy per bit.
- **Higher-order schemes pay for spectral efficiency.** At equal $E_b/N_0$, BER increases from BPSK/QPSK to 8PSK to 16-QAM, since more points are packed into the same energy budget.
- **Mapping matters.** The non-Gray QPSK curve sits above the Gray one, with the gap predicted by $\tfrac{3}{2}p - p^2$.
- **BFSK is 3 dB behind BPSK**, as theory predicts for orthogonal signalling.
- **Simulation tracks theory** within Monte Carlo error. The 8PSK and 16-QAM references are tight approximations, so small deviations at low SNR are expected.

## Known limitations

- **BER floor.** With $1.2 \times 10^6$ bits, the smallest measurable BER is about $8 \times 10^{-7}$. At high $E_b/N_0$ (BPSK and QPSK near 12 to 14 dB) the true BER is below this, so the simulation records zero errors and those points do not appear on the log plot. Increase `numBits` (keeping it a multiple of 12) to push the floor lower, at the cost of runtime and memory.
- **BFSK tone spacing.** The assignment specifies $f_i = (n_c+i)/T_b$, giving a spacing of $1/T_b$. Coherent orthogonality only requires $1/(2T_b)$, so this spacing is wider than the minimum and the BFSK signal occupies more bandwidth than necessary.

---

# Part 2: Rayleigh Flat Fading

*To be continued.*

**Scope:** BER of BPSK, QPSK and 16-QAM over a Rayleigh flat fading channel with perfect channel knowledge at the receiver.

System model:

$$
y_k = (h_r + j h_i)_k\, x_k + (n_c + j n_s)_k, \qquad h_r,\, h_i \sim \mathcal{N}\!\left(0, \tfrac{1}{2}\right)
$$

The receiver compensates for the channel gain, correlates, and makes a hard decision.

| Item | Status |
|------|--------|
| BPSK over Rayleigh | ☐ |
| QPSK over Rayleigh | ☐ |
| 16-QAM over Rayleigh | ☐ |
| Comparison vs. AWGN results from Part 1 | ☐ |

---

# Part 3: OFDM System

*To be continued.*

**Scope:** an OFDM link with the chain *coding → interleaver → mapper → 256-point IFFT → cyclic extension → channel → receiver*, evaluated over two channel models and with and without coding.

| Block | Specification |
|-------|---------------|
| Coding | none, or rate-1/5 repetition code (51 data bits per QPSK symbol plus one padding zero to fill 256) |
| Interleaver | $32 \times 16$ for QPSK, $32 \times 32$ for 16-QAM |
| Mapper | same as the single-carrier system |
| IFFT | 256 points |
| Channels | Rayleigh flat fading; frequency-selective fading with each sub-channel independently Rayleigh-faded in the frequency domain |

Required output: BER vs. SNR for BPSK and 16-QAM, for both channels, with and without coding (4 cases per modulation).

| Item | Status |
|------|--------|
| Transmitter chain | ☐ |
| Flat-fading receiver | ☐ |
| Frequency-selective receiver | ☐ |
| BER curves (4 cases per modulation) | ☐ |

---

## Getting Started

**Requirements**

- MATLAB R2019b or later
- Signal Processing Toolbox (`pwelch`, `hamming`)

**Run Part 1**

```matlab
cd part1
Part_1
```

The script runs end to end: it generates the data, simulates every scheme across the $E_b/N_0$ sweep, prints the BER table, and opens all nine figures.

**Tuning**

| Variable | Effect |
|----------|--------|
| `numBits` | more bits means a lower BER floor and smoother curves, but longer runtime |
| `EbNo_dB` | change the SNR sweep range or step |
| `numBits_PSD` | number of bits used for the BFSK Welch PSD |

## Repository Structure

Suggested layout:

```
ber-analysis-awgn-fading-ofdm/
├── README.md
├── part1/                  # Single carrier over AWGN
│   └── Part_1.m
├── part2/                  # Rayleigh flat fading
├── part3/                  # OFDM
├── docs/                   # Assignment PDFs, reports
└── results/
    ├── part1/              # Exported figures
    ├── part2/
    └── part3/
```

## Team

| Part | Owner |
|------|-------|
| 1: Single-carrier / AWGN | Mohamed Aboualy ([GitHub](https://github.com/mo-aboualy)) |
| 2: Rayleigh flat fading | *Teammate name* |
| 3: OFDM | *Teammate name* |

## Acknowledgements

Course project for Digital Communications, Faculty of Engineering, Cairo University.
