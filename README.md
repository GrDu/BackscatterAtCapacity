# BackscatterLoadModulationAtCapacity
 Matlab code associated with the publication:
 "Load Modulation for Backscatter Communication: Channel Capacity and Capacity-Approaching Finite Constellations"
 arXiv preprint arXiv:2207.08100, July 2022.  
 Available online at (open access): https://arxiv.org/abs/2207.08100
 
The numbers in the filenames relate to the figure numbers in the paper. For example, the script eval03_*.m generates Figure 3. The generated Matlab figures sometimes hold more content and have slightly different style (b/c the paper figures use tikz).

A secondary merit of this repository are the provided functions that return the data rate over certain additive white Gaussian noise (AWGN) channels. These can find use in all kinds of communications analyses. The specific functionality comprises:
- Compute the achievable information rate over an AWGN channel (real- or complex-valued) with given noise variance, for a given finite transmit symbol alphabet (i.e. a given constellation) with given symbol probabilities.
- Compute the achievable information rate  over a complex-valued AWGN channel for a transmit signaling with given discrete amplitudes and uniform independent phase (DAUIP).
- Compute the channel capacity of a peak-power-limited AWGN channel (real- or complex-valued).
