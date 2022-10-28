# BackscatterLoadModulationAtCapacity
 Matlab code associated with the publication:
 "Load Modulation for Backscatter Communication: Channel Capacity and Capacity-Approaching Finite Constellations"
 arXiv preprint arXiv:2207.08100, July 2022.  
 Available online at (open access): https://arxiv.org/abs/2207.08100
 
The numbers in the filenames relate to the figure numbers in the paper. For example, the script eval03_*.m generates Figure 3. The generated Matlab figures sometimes hold more content and have slightly different style (b/c the paper figures use tikz).

Secondary scientific utility of this repository stems from the provided functions concerning the information theory of additive white Gaussian noise (AWGN) channels in real- and complex-valued (i.e. quadrature) guise. These can be useful in all kinds of communications analyses. The specific functionality comprises:
- Compute the achievable information rate of a given finite symbol alphabet, with given symbol probabilities, over an AWGN channel (real- or complex-valued) with given noise variance.
- Compute the channel capacity of a peak-power-limited AWGN channel (real- or complex-valued).
- Compute the achievable information rate with discrete-amplitude-and-uniform-independent-phase (DAUIP) transmit signaling over a complex-valued AWGN channel.
