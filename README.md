## Sparsity-Assisted Signal Denoising and Pattern Recognition
Signal denoising and pattern recognition of time-series data are widely used in many scientific fields, including physics, engineering, medicine, economics, acoustics, biology, and psychology. For example, specific signal patterns in the electroencephalogram (EEG) data, which are useful in clinical diagnosis and cognitive neuroscience, are challenging to detect and distinguish from artifacts. In this work, we address the problem of signal denoising and pattern recognition in processing batch-mode time-series data by combining linear time-invariant filters, orthogonal multiresolution representations, and sparsity-based methods. For example, in the figure below, we decompose the original signal as the sum of a low-frequency signal, an oscillatory signal, and a discontinuous signal using our novel filter designs and proposed signal model.

<p align="center">
  <img width="500" height="450" src="https://github.com/prateekgv/sasdpr/blob/master/images/sasdpr.png">
</p>

### Main Contributions

#### Novel Zero-Phase IIR Filters as Matrices
* The zero-phase filters designed as matrices are stable and not sparse.
* The orders of the filter depend on the positive definiteness condition of the reachability and observability Gramians.
* Filter response types include stable low-pass, high-pass, and band-pass filters.

#### Applications
* We develop a new signal model called sparsity-assisted signal denoising (SASD) by combining our proposed filter designs with the existing signal model. Because the zero-phase filters in the SASD signal model are stable, they demonstrate consistent results on changing the orders of the filter.
* We propose and derive a new signal model called sparsity-assisted pattern recognition (SAPR). In SAPR, we combine LTI band-pass filters and sparsity-based methods with orthogonalmultiresolution representations, such as wavelets, to detect specific patterns in the input signal.
* We combine the signal denoising and pattern recognition tasks, and derive a new signal model called the sparsity-assisted signal denoising and pattern recognition (SASDPR).

#### Cite
If you find our work interesting and useful, please cite the following publication:
> G. V. Prateek, Y-E. Ju, and A. Nehorai, “Sparsity-assited signal denoising and pattern recognition in time-series data,” to appear in _Circuits, Systems, and Signal Processing_.
