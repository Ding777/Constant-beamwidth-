
#  Circular Microphone Array MVDR Beamforming with Diffuse Noise Field and Sidelobe Control

## 1. Overview

This MATLAB script implements **frequency-domain beamforming** for a **uniform circular microphone array (UCA)** under a **diffuse noise field**.

The code performs:

* Modeling of **diffuse noise spatial coherence**
* **MVDR (Minimum Variance Distortionless Response)** beamformer design
* Beam pattern visualization over **frequency × azimuth**
* **Sidelobe-constrained optimization** using **CVX**
* Comparison between:

  * Standard MVDR beamformer
  * Optimized sidelobe-controlled beamformer
  * Re-weighted MVDR solution

This script is mainly for **array signal processing research**, **spatial filtering**, and **beam pattern analysis**.

---

## 2. Array and Signal Parameters

```matlab
f01 = 4000;        % Signal center frequency (Hz)
B = 8000;          % Bandwidth (Hz)
Nfft = 512;        % FFT size
```

* Frequency bins are defined from `0 → Fs/2`
* `f_range` maps FFT bins to **physical frequencies**

```matlab
M = 16;            % Number of microphones
c = 340;           % Speed of sound (m/s)
R = 0.025;         % Radius of circular array (meters)
```

The array is a **16-element Uniform Circular Array (UCA)**.

Each microphone angle:

```matlab
gamma = 2*pi*[0:M-1]'/M;
```

---

## 3. Source Direction (DOA)

```matlab
seta0 = 60*pi/180;     % Elevation angle
fai0  = -150*pi/180;  % Azimuth angle
```

* Defines the **desired signal direction**
* Used to compute the **steering vector**

---

## 4. Diffuse Noise Field Modeling

### 4.1 Spatial Coherence Matrix

```matlab
Fvv(:,i,j) = sin(2*pi*f0*dij/c)./(2*pi*f0*dij/c);
```

* Models **diffuse noise spatial coherence**
* Based on **sinc function**
* `dij` = distance between microphone `i` and `j`

Properties:

* Diagonal = 1 (self-coherence)
* Off-diagonal = frequency-dependent coherence

This is a **standard diffuse noise model**.

---

### 4.2 Diagonal Loading

```matlab
Fvv2(n,:,:) = inv(Fvv(n,:,:) + eye(M)*0.1);
```

* Adds diagonal loading for **numerical stability**
* Improves robustness of MVDR
* Prevents matrix singularity

---

## 5. MVDR Beamformer Design

### 5.1 Steering Vector (UCA)

```matlab
a0 = exp(-1j*2*pi*R*freq/c*sin(seta0)*cos(fai0-gamma));
```

This is the **exact UCA steering vector**, accounting for:

* Array radius
* Elevation angle
* Azimuth angle

---

### 5.2 MVDR Weight Computation

```matlab
Wc = (Rxx * a0) / (a0' * Rxx * a0);
```

This enforces:

* **Distortionless response** in desired direction
* **Minimum output power**

This is the **classical MVDR beamformer**.

---

## 6. Beam Pattern Scanning

### 6.1 Azimuth Scan

```matlab
scan_angle = [-180:1:180]/180*pi;
```

For each angle:

* Compute steering vector
* Evaluate beam response

```matlab
p = 20*log10(abs(a_all' * conj(Wopt)) / max(...));
```

This produces:

* **Azimuthal beam pattern (dB)**
* Normalized to peak = 0 dB

---

### 6.2 Frequency-Angle Beam Pattern

```matlab
surf(angle, frequency, beam_response)
```

This visualizes:

* Beam width vs frequency
* Frequency-dependent sidelobes
* Array resolution behavior

---

## 7. Mainlobe and Sidelobe Region Definition

At a reference frequency:

```matlab
angle_ml = [-180:-110]/180*pi;   % Main lobe
angle_sl = remaining angles;    % Side lobes
```

Purpose:

* Separate **desired region** and **interference region**
* Required for sidelobe optimization

---

## 8. Sidelobe-Constrained Optimization (CVX)

### 8.1 Optimization Goal

```matlab
minimize || p_mainlobe - desired_response ||_2
```

### 8.2 Constraints

```matlab
abs(p_sidelobe) <= max_sl
w' * a0 == 1
||w||_2 <= M
```

Interpretation:

* Keep mainlobe shape
* Suppress sidelobes
* Preserve unity gain in look direction
* Control weight norm

This is a **convex beamforming optimization**.

---

## 9. Post-Optimization MVDR Reweighting

```matlab
Wopt1 = (Rxx * w_msl) / (w_msl' * Rxx * w_msl);
```

* Uses optimized weights as steering constraint
* Combines **MVDR robustness** with **sidelobe control**

---

## 10. Visualization Outputs

### Figures Produced

| Figure | Description                                 |
| ------ | ------------------------------------------- |
| Fig.1  | MVDR beam patterns (multiple frequencies)   |
| Fig.3  | Optimized sidelobe-controlled beam patterns |
| Fig.4  | Optimized frequency-angle surface           |
| Fig.5  | Original MVDR frequency-angle surface       |

---

## 11. Key Techniques Used

* Uniform Circular Array (UCA)
* Diffuse noise spatial coherence
* MVDR beamforming
* Diagonal loading
* Frequency-domain beam scanning
* Convex optimization (CVX)
* Sidelobe suppression
* Frequency-dependent beam analysis

---

## 12. Intended Applications

This code is suitable for:

* Microphone array research
* Speech enhancement
* Spatial audio processing
* DOA-aware beamforming
* Academic demonstrations of MVDR + optimization

---

## 13. Dependencies

* MATLAB
* **CVX Toolbox** (required for sidelobe optimization)


