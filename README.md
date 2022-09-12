# BackscatterLoadModulationAtCapacity
 Matlab code associated with the publication:
 "Load Modulation for Backscatter Communication: Channel Capacity and Capacity-Approaching Finite Constellations"
 arXiv preprint arXiv:2207.08100, July 2022.  
 Available online at (open access): https://arxiv.org/abs/2207.08100
 
The numbers in the filenames relate to the figure number in the paper. For example, the script eval03_*.m generates Figure 3. The figures in the paper have slightly different style (tikz) and more tidy content selection.

Apart from backscatter load modulation, the repository contains useful information-theoretic functions concerning additive white Gaussian noise (AWGN) channels in real- and complex-valued (i.e. quadrature) guise. In particular:
- Computing the achievable information rate with a finite symbol alphabet over an AWGN channel (real- or complex-valued).
- Computing the achievable information rate with discrete-amplitude-and-uniform-independent-phase (DAUIP) transmit signaling over a complex-valued AWGN channel.
- Computing the channel capacity of a peak-power-limited AWGN channel (real- or complex-valued).
