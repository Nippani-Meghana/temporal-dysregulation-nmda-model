# Modeling Excitatory–Inhibitory Dynamics in Schizophrenia with a Minimal Morris–Lecar Microcircuit

## Abstract

This repository implements a biophysically minimal two-neuron microcircuit (one excitatory, one inhibitory) using the Morris–Lecar formalism to examine how NMDA receptor hypofunction disrupts excitation–inhibition (E/I) balance and temporal stability. The model couples an excitatory (E) unit to an inhibitory (I) unit via NMDA (E→I) and GABA\_A (I→E) synapses with distinct kinetics (τ\_NMDA ≫ τ\_GABA). We compare a "Normal" condition against an "NMDA hypofunction" condition and quantify timing stability using inter-spike interval (ISI) statistics, coefficient of variation (CV), firing rates (FR), and the network-level E/I rate ratio. The minimal circuit reveals a characteristic signature: inhibitory irregularity increases while excitatory regularity is preserved, tilting the E/I balance toward excitation and degrading temporal precision.

---

## Important Cautionary Note

This project is an **independent research exercise** in computational modeling. It is not affiliated with any university, laboratory, or clinical institution. The values, parameters, and results produced by this code are **not calibrated to clinical or biological data** and should not be interpreted as reflecting actual patient physiology or diagnostic outcomes. The model is intended **solely for educational and exploratory purposes** in computational neuroscience.

---

## Key Contributions

* **Minimal, transparent microcircuit**: Two coupled Morris–Lecar units with slow NMDA and fast GABA synapses.
* **Mechanistic manipulation**: Isolated reduction of NMDA gain captures schizophrenia-like inhibitory weakening and timing irregularity.
* **Quantitative readouts**: ISI distributions, CV(ISI), firing rates, and E/I ratio comparisons (Normal vs Schizo).
* **Reproducible figures**: Spike traces and ISI histograms, plus summary bar plots for CV, FR, and E/I ratio.

---

## Model Overview

* **Single-compartment Morris–Lecar neurons** with standard voltage-dependent calcium activation \$m\_\infty(V)\$ and potassium recovery \$w\$ dynamics.
* **Synaptic coupling**:

  * E→I via **NMDA**: slow gating \$s\_E\$ with logistic activation and \$\tau\_{NMDA} \approx 100 , ms\$ (model units), reversal \$E\_{exc} = 0\$.
  * I→E via **GABA\_A**: fast gating \$s\_I\$ with \$\tau\_{GABA} \approx 10 , ms\$, reversal \$E\_{inh} \in \[-0.8,-0.7]\$.
* **Integration**: Explicit Euler with time step `dt = 0.01` (model units). A detector refractory of \~8 ms is used to avoid double-counting threshold crossings.
* **Noise**: Optional additive Gaussian noise to membrane voltage per step to emulate stochastic drive; reduced in Normal, slightly higher in Schizo.

> **Caution on units and interpretation:** The membrane potential and time scales here are expressed in **model units** (e.g., \$V\_{Ca} = 1, V\_K = -0.7\$). Axis labels therefore use **a.u.** for voltage. If biological millivolts/milliseconds are desired, apply an explicit rescaling. In the present configuration, firing rates and ISIs reflect model dynamics rather than directly calibrated cortical values; conclusions are drawn **comparatively** (Normal vs Schizo) rather than absolutely.

---

## Conditions Compared

* **Normal** (baseline): Higher NMDA conductance, balanced tonic drive to E and I, moderate GABA.
* **Schizo (NMDA hypofunction)**: Reduced NMDA conductance (E→I), slightly altered tonic inputs and inhibitory coupling to reflect hypofunction and increased background variability.

Parameters are defined in the script blocks `normal_parameters` and `Schizophrenia_parameters`. For clean attribution of effects, we recommend varying **one parameter at a time** (e.g., only `g_nmda`) when performing mechanistic analyses.

---

## Outputs and Diagnostics

1. **Spike Traces**

   * Overlaid E and I voltage trajectories with spike times (for Normal and Schizo).
2. **ISI Distributions**

   * Histograms of Excitatory and Inhibitory ISIs (Normal vs Schizo, shared bin edges).
3. **Summary Metrics**

   * **CV(ISI)**: \$CV = \sigma(ISI)/\mu(ISI)\$ — irregularity index.
   * **Firing Rate (FR)**: \$FR = \text{spikes}/\text{duration}\$.
   * **E/I Ratio**: \$FR\_E / FR\_I\$ (network-level balance).
4. **Comparison Plots**

   * Bar plots for CV(ISI), FR, and E/I ratio across conditions.

**Interpretation template** (observed in typical runs):

* Excitatory ISIs remain narrowly concentrated (clock-like), CV low.
* Inhibitory ISIs broaden under NMDA hypofunction, CV increases, and occasional cycle skips appear.
* FR of I may drop slightly relative to E, pushing **E/I ratio > 1** in Schizo.

---

## How to Run

1. Install dependencies:

   * Python 3.9+
   * `numpy`, `matplotlib`
2. Execute the script (e.g., `python main.py`). The script will:

   * Simulate **Normal** and **Schizo** conditions,
   * Display spike traces and ISI histograms,
   * Compute and plot **CV(ISI)**, **FR**, and **E/I ratio**.

If adapting to notebooks, ensure plots render inline and that random seeds are controlled for exact reproducibility when needed.

---

## Reproducibility Notes

* **Randomness**: The noise term is stochastic. Set a fixed seed (e.g., `np.random.seed(...)`) for reproducible spike times and derived statistics.
* **Detector threshold vs synaptic thresholds**: Spike detection threshold is independent of synaptic logistic centers; coupling them can bias counts. Keep them distinct for clean analyses.
* **Parameter sweeps**: To attribute causality, sweep one parameter (e.g., `g_nmda`) across a range while holding others constant and report FR, CV, and E/I ratio as functions of that parameter.

---

## Limitations and Scope

* **Minimality**: Two units cannot capture population phenomena (e.g., network oscillations, heterogeneity). Results emphasize **directional** effects rather than absolute magnitudes.
* **Biophysical calibration**: Without explicit unit mapping, voltages are in **a.u.** and time constants reflect model scaling. If biological realism is required, calibrate to literature benchmarks and adjust \$C\$, conductances, and `rate_const` to achieve cortical firing ranges (E ≈ 5–30 Hz; I ≈ 20–80 Hz) before interpretation.
* **Noise modeling**: Additive voltage noise is a proxy for synaptic/channel variability; in larger models, input current noise or Poisson synaptic bombardment may be preferable.

---


## Ethical and Terminology Note

“Schizo” in code/variable names is shorthand for an NMDA-hypofunction condition used in computational modeling. 
