# BER Analysis of Digital Modulation Schemes over AWGN, Rayleigh Flat-Fading, and Frequency-Selective Channels with OFDM

**A from-scratch MATLAB simulation study of digital modulation, from single-carrier AWGN links to Rayleigh fading and a full OFDM transceiver.**

Digital Communications · 3rd Year · Electronics & Communications Dept. · Faculty of Engineering, Cairo University

![MATLAB](https://img.shields.io/badge/MATLAB-R2019b%2B-orange)
![Toolbox](https://img.shields.io/badge/Part%201%20requires-Signal%20Processing%20Toolbox-blue)
![Status](https://img.shields.io/badge/status-Parts%201%E2%80%932%20done%20%7C%20Part%203%20in%20progress-yellow)

---

## Overview

This repository is a three-part simulation project that builds a digital communication link step by step. Every block (mapper, channel, demapper, BER counter) is written by hand in plain MATLAB with no `comm.*` toolbox objects, and every simulated curve is checked against a closed-form theoretical BER.

| Part | Topic | Channel | Status |
|------|-------|---------|--------|
| **1** | Single-carrier system: BPSK, QPSK (Gray and non-Gray), 8PSK, 16-QAM, BFSK, plus BFSK PSD | AWGN | Done |
| **2** | BER of BPSK, QPSK and 16-QAM with channel inversion at the receiver | Rayleigh flat fading, channel known at the receiver | Done |
| **3** | OFDM transceiver: coding, interleaver, mapper, 256-point IFFT, cyclic extension | Rayleigh flat fading and frequency-selective fading | In progress |

> **Deadline:** all parts due 14 Dec 2026, 11:59 PM.

---

## Table of Contents

1. [Part 1: Single-Carrier System](#part-1-single-carrier-system)
2. [Part 2: Rayleigh Flat Fading](#part-2-rayleigh-flat-fading)
3. [Part 3: OFDM System](#part-3-ofdm-system) *(to be continued)*
4. [Getting Started](#getting-started)
5. [Repository Structure](#repository-structure)

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
| 16-QAM, Gray | 4 | $\lbrace 00, 01, 11, 10 \rbrace \to \lbrace -3, -1, +1, +3 \rbrace$ per axis | 10 | 2.5 | per-axis thresholds at $0$ and $\pm 2$ (ML) |
| BFSK | 1 | orthonormal basis, $0 \to [\sqrt{E_b}, 0],\ 1 \to [0, \sqrt{E_b}]$ | 1 | 1 | correlator comparison $r_2 > r_1$ |

### Noise scaling

For every scheme the noise power is set from the target $E_b/N_0$ using that scheme's own $E_b$:

$$
N_0 = \frac{E_b}{E_b/N_0}, \qquad
n \sim \mathcal{N}\left(0, \frac{N_0}{2}\right) \text{ per real dimension}
$$

This keeps the comparison across modulation orders fair, since $E_b$, not $E_s$, is the common reference.

## Theoretical references

Simulated points are overlaid on these expressions, with $Q(x) = \frac{1}{2}\mathrm{erfc}\left(x/\sqrt{2}\right)$:

| Scheme | Theoretical BER | Type |
|--------|-----------------|------|
| BPSK | $\frac{1}{2}\mathrm{erfc}\left(\sqrt{E_b/N_0}\right)$ | exact |
| QPSK, Gray | $\frac{1}{2}\mathrm{erfc}\left(\sqrt{E_b/N_0}\right)$ | exact |
| QPSK, non-Gray | $\frac{3}{2}p - p^2,\quad p = \frac{1}{2}\mathrm{erfc}\left(\sqrt{E_b/N_0}\right)$ | exact |
| 8PSK, Gray | $\frac{1}{3}\mathrm{erfc}\left(\sqrt{3E_b/N_0}\,\sin\frac{\pi}{8}\right)$ | tight approximation |
| 16-QAM, Gray | $\frac{3}{8}\mathrm{erfc}\left(\sqrt{\frac{2}{5}\,E_b/N_0}\right)$ | tight approximation |
| BFSK, coherent orthogonal | $\frac{1}{2}\mathrm{erfc}\left(\sqrt{\frac{E_b}{2N_0}}\right)$ | exact |

**Why non-Gray QPSK is worse.** In Fig. 3, a decision error on the $I$ axis flips both bits (because $b_1$ depends on $I$), while an error on $Q$ alone flips one bit, and errors on both axes flip one bit net. Averaging over the two bits gives $\frac{3}{2}p - p^2$, which is about $1.5\times$ the Gray-coded BER at high SNR. Gray mapping ensures nearest-neighbour symbol errors cost exactly one bit error.

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

<!-- After running the script, export the figures into part1-awgn/figures/ and uncomment:

![All schemes](part1-awgn/figures/fig5_all_schemes.png)
![QPSK Gray vs non-Gray](part1-awgn/figures/fig6_qpsk_gray_vs_nongray.png)
![BFSK PSD](part1-awgn/figures/fig8_bfsk_psd.png)
-->

## Key observations

- **BPSK and Gray QPSK overlap exactly** on a BER vs. $E_b/N_0$ plot: QPSK is two independent BPSK streams on orthogonal axes with the same energy per bit.
- **Higher-order schemes pay for spectral efficiency.** At equal $E_b/N_0$, BER increases from BPSK/QPSK to 8PSK to 16-QAM, since more points are packed into the same energy budget.
- **Mapping matters.** The non-Gray QPSK curve sits above the Gray one, with the gap predicted by $\frac{3}{2}p - p^2$.
- **BFSK is 3 dB behind BPSK**, as theory predicts for orthogonal signalling.
- **Simulation tracks theory** within Monte Carlo error. The 8PSK and 16-QAM references are tight approximations, so small deviations at low SNR are expected.

## Known limitations

- **BER floor.** With $1.2 \times 10^6$ bits, the smallest measurable BER is about $8 \times 10^{-7}$. At high $E_b/N_0$ (BPSK and QPSK near 12 to 14 dB) the true BER is below this, so the simulation records zero errors and those points do not appear on the log plot. Increase `numBits` (keeping it a multiple of 12) to push the floor lower, at the cost of runtime and memory.
- **BFSK tone spacing.** The assignment specifies $f_i = (n_c+i)/T_b$, giving a spacing of $1/T_b$. Coherent orthogonality only requires $1/(2T_b)$, so this spacing is wider than the minimum and the BFSK signal occupies more bandwidth than necessary.

---

# Part 2: Rayleigh Flat Fading

Part 2 replaces the AWGN-only channel of Part 1 with a **Rayleigh flat-fading channel** and measures the BER of **BPSK, QPSK and 16-QAM** when the receiver has **perfect knowledge of the channel** and compensates for it by channel inversion.

```mermaid
flowchart LR
    A[Random bits] --> B[Mapper]
    B --> C["Multiply by h (Rayleigh)"]
    C --> D["Add n (AWGN)"]
    D --> E[Channel inversion]
    E --> F[Hard-decision demapper]
    F --> G[BER calculator]
    A -.-> G
```

## System model

Each transmitted symbol $x_k$ is scaled and rotated by an independent complex channel coefficient $h_k$ and corrupted by complex noise $n_k$:

$$
y_k = (h_r + j h_i)_k \, x_k + (n_c + j n_s)_k
$$

| Quantity | Distribution | Meaning |
|----------|--------------|---------|
| $h_r,\ h_i$ | $\mathcal{N}\left(0, \frac{1}{2}\right)$, independent | $h = A_h e^{j\theta_h}$ with Rayleigh amplitude $A_h$ and uniform phase $\theta_h \in [0, 2\pi)$ |
| $n_c,\ n_s$ | $\mathcal{N}\left(0, \frac{N_0}{2}\right)$, independent | complex AWGN |
| $\mathbb{E}\left[\lvert h \rvert^2\right]$ | $1$ | average channel power gain is normalized, so $E_b/N_0$ is the **average** SNR per bit |

**Receiver (channel known).** The received sample is divided by the channel coefficient and a hard decision with threshold 0 is taken:

$$
\tilde{y}_k = \frac{y_k}{h_k} = x_k + \frac{n_k}{h_k}
$$

After inversion the noise on each real axis has variance $\frac{N_0}{2\lvert h_k \rvert^2}$, so a deep fade (small $\lvert h_k \rvert$) amplifies the noise. This noise enhancement, not the average noise level, is what dominates the error rate over fading channels.

**Noise scaling** follows Part 1: $N_0 = E_b / (E_b/N_0)$ using each scheme's own $E_b$.

## What is implemented

| Scheme | Bits/symbol | Mapping | $E_s$ | $E_b$ | Hard decision on $\tilde{y}$ |
|--------|:-----------:|---------|:-----:|:-----:|------------------------------|
| BPSK | 1 | $0 \to -\sqrt{E_b},\ 1 \to +\sqrt{E_b}$ | 1 | 1 | $\hat{b} = \mathbb{1}\left[\mathrm{Re}\tilde{y} > 0\right]$ |
| QPSK | 2 | $b_0 \to I,\ b_1 \to Q$, each $\lbrace 0, 1 \rbrace \to \lbrace -1, +1 \rbrace$ | 2 | 1 | independent sign decisions on $\mathrm{Re}\tilde{y}$ and $\mathrm{Im}\tilde{y}$ |
| 16-QAM, Gray | 4 | $b_0 b_1 \to I,\ b_2 b_3 \to Q$, each $\lbrace 00, 01, 11, 10 \rbrace \to \lbrace -3, -1, +1, +3 \rbrace$ | 10 | 2.5 | per axis: sign gives $b_0$ (or $b_2$), $\lvert\cdot\rvert < 2$ gives $b_1$ (or $b_3$) |

QPSK and 16-QAM follow the constellations given in the assignment. The 16-QAM thresholds are at $0$ and $\pm 2$, the midpoints between adjacent levels.

## Simulation procedure

The code follows steps (a) to (h) of the assignment for every modulation:

| Step | Description | Implementation |
|:----:|-------------|----------------|
| a | random bit stream | `randi([0 1], 1, N_bits)`, shared by all three schemes |
| b | symbol mapping | polar NRZ (BPSK), I/Q pairs (QPSK), Gray levels (16-QAM) |
| c | channel and noise vectors | $h_r, h_i \sim \mathcal{N}(0, 1/2)$ and $n_c, n_s \sim \mathcal{N}(0, N_0/2)$, new draws for every symbol and every SNR point |
| d | received symbols | `y = h .* x + n` |
| e | channel compensation and decision | `y ./ h`, then thresholding (see the table above) |
| f | BER | `sum(b ~= b_hat) / N_bits` |
| g | sweep the noise level | loop over $E_b/N_0$ |
| h | plot BER vs. SNR | `semilogy`, one figure per scheme plus comparisons |

The channel, noise generation and inversion steps are shared by the three schemes, so they are written once in a helper function (`pass_through_rayleigh_channel`).

| Parameter | Value | Note |
|-----------|-------|------|
| Bits per $E_b/N_0$ point | $10^6$ | divisible by 1, 2 and 4 bits/symbol with no padding |
| $E_b/N_0$ range | $0$ to $30$ dB, 2 dB steps | wider than Part 1 because fading needs much more SNR |
| Fading | independent $h_k$ for every symbol | i.i.d. flat fading (ideal interleaving); the average BER matches the closed forms below |
| Random seed | not fixed | curves vary slightly between runs, within Monte Carlo error |

## Theoretical references

Averaging the AWGN error probability over the exponential distribution of $\lvert h \rvert^2$ uses the identity

$$
\mathbb{E}\left[Q\left(\sqrt{c\,\lvert h \rvert^2}\right)\right] = \frac{1}{2}\left(1 - \sqrt{\frac{c}{2 + c}}\right), \qquad \lvert h \rvert^2 \sim \text{Exp}(1)
$$

Writing $\gamma = E_b/N_0$ for the average SNR per bit:

| Scheme | BER over Rayleigh fading (coherent, channel known) | High-SNR behaviour |
|--------|-----------------------------------------------------|--------------------|
| BPSK and QPSK | $\displaystyle \frac{1}{2}\left(1 - \sqrt{\frac{\gamma}{1+\gamma}}\right)$ | $\approx \dfrac{1}{4\gamma}$ |
| 16-QAM, Gray | $\displaystyle \frac{1}{8}\left[4 - 3\sqrt{\frac{0.4\gamma}{1+0.4\gamma}} - 2\sqrt{\frac{3.6\gamma}{1+3.6\gamma}} + \sqrt{\frac{10\gamma}{1+10\gamma}}\right]$ | $\approx \dfrac{0.497}{\gamma}$ |

<details>
<summary>Derivation of the 16-QAM expression</summary>

With Gray mapping, 16-QAM is two independent 4-PAM streams with levels $\pm 1, \pm 3$ on $I$ and $Q$. For a fixed $\lvert h \rvert^2$ the equalized noise per axis has variance $\frac{N_0}{2\lvert h \rvert^2}$. Let $a = \sqrt{2\lvert h \rvert^2 / N_0}$ and recall $N_0 = 2.5/\gamma$, so $a^2 = 0.8\,\gamma\,\lvert h \rvert^2$.

- Sign bit ($b_0$): error probability $\frac{1}{2}\left[Q(a) + Q(3a)\right]$.
- Magnitude bit ($b_1$): error probability $\frac{1}{2}\left[2Q(a) + Q(3a) - Q(5a)\right]$.

Averaging the two bits gives the conditional BER

$$
P_b(\lvert h \rvert^2) = \frac{1}{4}\left[3Q(a) + 2Q(3a) - Q(5a)\right]
$$

Applying the identity above with $c = 0.8\gamma k^2$ for $k = 1, 3, 5$ yields the closed form in the table. At high SNR, $\gamma \to \infty$ expands each square root as $1 - \frac{1}{2x}$ and gives $\approx 0.497/\gamma$.

</details>

The AWGN reference curves overlaid in the comparison figures are the Part 1 expressions: $\frac{1}{2}\mathrm{erfc}\left(\sqrt{\gamma}\right)$ for BPSK and QPSK, and $\frac{3}{8}\mathrm{erfc}\left(\sqrt{\frac{2}{5}\gamma}\right)$ for 16-QAM.

**Theoretical Rayleigh BER at selected points**

| $E_b/N_0$ | BPSK / QPSK | 16-QAM |
|:---------:|:-----------:|:------:|
| 0 dB | $1.46 \times 10^{-1}$ | $1.98 \times 10^{-1}$ |
| 10 dB | $2.33 \times 10^{-2}$ | $4.24 \times 10^{-2}$ |
| 20 dB | $2.48 \times 10^{-3}$ | $4.89 \times 10^{-3}$ |
| 30 dB | $2.50 \times 10^{-4}$ | $4.96 \times 10^{-4}$ |

## Results

**Per-scheme figures:** simulated Rayleigh BER (left) and the same curve against the AWGN reference (right).

| Scheme | Rayleigh BER | Rayleigh vs. AWGN |
|--------|--------------|-------------------|
| BPSK | ![BPSK over Rayleigh](part2-fading/results/BER%20Performance%20of%20BPSK%20over%20Rayleigh%20Flat%20Fading%20Channel.jpg) | ![BPSK Rayleigh vs AWGN](part2-fading/results/BPSK%20BER%20Comparison%20Rayleigh%20Flat%20Fading%20vs.%20AWGN%20Channel.jpg) |
| QPSK | ![QPSK over Rayleigh](part2-fading/results/BER%20Performance%20of%20QPSK%20over%20Rayleigh%20Flat%20Fading%20Channel.jpg) | ![QPSK Rayleigh vs AWGN](part2-fading/results/QPSK%20BER%20Comparison%20Rayleigh%20Flat%20Fading%20vs.%20AWGN%20Channel.jpg) |
| 16-QAM | ![16-QAM over Rayleigh](part2-fading/results/BER%20Performance%20of%2016-QAM%20over%20Rayleigh%20Flat%20Fading%20Channel.jpg) | ![16-QAM Rayleigh vs AWGN](part2-fading/results/16-QAM%20BER%20Comparison%20Rayleigh%20Flat%20Fading%20vs.%20AWGN%20Channel.jpg) |

**Overall comparison over the Rayleigh channel**

![BPSK, QPSK and 16-QAM over Rayleigh](part2-fading/results/BER%20Comparison%20BPSK%2C%20QPSK%2C%20and%2016-QAM%20over%20Rayleigh%20Channel.jpg)

Running the script produces seven figures:

| # | Figure |
|---|--------|
| 1 | BPSK over Rayleigh fading |
| 2 | BPSK: Rayleigh vs. AWGN |
| 3 | QPSK over Rayleigh fading |
| 4 | QPSK: Rayleigh vs. AWGN |
| 5 | 16-QAM over Rayleigh fading |
| 6 | 16-QAM: Rayleigh vs. AWGN |
| 7 | BPSK, QPSK and 16-QAM over Rayleigh fading |

## Key observations

- **Fading is far more costly than noise alone.** To reach a BER of $10^{-3}$, BPSK and QPSK need about 24.0 dB over Rayleigh fading versus 6.8 dB over AWGN, a penalty of roughly 17 dB.
- **The curves are straight lines with slope $-1$ on the log-log scale.** BER falls only 10 times per 10 dB ($\approx 1/(4\gamma)$), compared with the exponential decay over AWGN. Errors are caused mostly by rare deep fades where $\lvert h \rvert$ is small, not by the average noise level. This is the diversity-order-1 behaviour of a single-antenna flat-fading link.
- **BPSK and QPSK still coincide** when plotted against $E_b/N_0$: each QPSK axis is an independent BPSK decision at the same energy per bit. The small difference visible near 30 dB is Monte Carlo error (about 250 errors per point at $10^6$ bits).
- **16-QAM needs about 3 dB more than BPSK/QPSK.** For $10^{-3}$ it needs about 26.9 dB, versus 24.0 dB for BPSK. The two curves are parallel, so the gap is a constant SNR offset rather than a change in slope. Over AWGN the gap is slightly larger, about 3.7 dB (6.8 dB vs. 10.5 dB, using the 16-QAM approximation).
- **Simulation tracks theory.** Simulated BER agrees with the closed forms above within Monte Carlo error across the whole sweep.

## Known limitations

- **Perfect channel knowledge.** The receiver divides by the exact $h_k$. A practical receiver would need channel estimation, which adds estimation error and pilot overhead.
- **BER floor and statistical accuracy.** With $10^6$ bits the smallest measurable BER is $10^{-6}$. At 30 dB the BER is about $2.5 \times 10^{-4}$, i.e. roughly 250 errors, so the points near the end of the sweep carry about $\pm 6\%$ (one standard deviation) statistical uncertainty. Increase `N_bits` (keeping it divisible by 4) to tighten them.
- **Memoryless fading.** $h_k$ is redrawn for every symbol (fast fading with ideal interleaving). Real channels are correlated in time, which changes burst behaviour but not the average BER for coherent detection with known $h$.

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

- MATLAB R2019b or later (local functions inside scripts are used in Part 2)
- Signal Processing Toolbox for Part 1 (`pwelch`, `hamming`); Part 2 needs no toolbox

**Run Part 1**

```matlab
cd part1-awgn
Part_1
```

The script runs end to end: it generates the data, simulates every scheme across the $E_b/N_0$ sweep, prints the BER table, and opens all nine figures.

**Run Part 2**

```matlab
cd part2-fading
part2_ber_analysis
```

The script simulates BPSK, QPSK and 16-QAM over the Rayleigh channel and opens all seven figures. Export them into `results/` to refresh the images shown above.

**Tuning**

| Variable | Script | Effect |
|----------|--------|--------|
| `numBits` | Part 1 | more bits means a lower BER floor and smoother curves, but longer runtime |
| `EbNo_dB` | Part 1 | change the SNR sweep range or step |
| `numBits_PSD` | Part 1 | number of bits used for the BFSK Welch PSD |
| `N_bits` | Part 2 | number of bits per SNR point (keep it divisible by 4) |
| `SNR_dB` | Part 2 | change the SNR sweep range or step |

## Repository Structure

```
ber-analysis-awgn-fading-ofdm/
├── README.md
├── part1-awgn/                    # Single carrier over AWGN
│   ├── Part_1.m
│   └── figures/                   # 9 result figures (.png)
├── part2-fading/                  # Rayleigh flat fading
│   ├── part2_ber_analysis.m
│   └── results/                   # 7 result figures (.jpg)
├── part3-ofdm/                    # OFDM
│   └── ofdm_project_complete.m
└── docs/                          # Assignment reports
    ├── Part1_Report.docx
    ├── Part1_Report.pdf
    ├── Part3_Report_V3.docx
    └── Part3_Report_V3.pdf
```

## Acknowledgements

Course project for Digital Communications, Faculty of Engineering, Cairo University.
