# Title:	Matlab code for "Computational Adaptive Optics for Fluorescence Microscopy via Sparse Blind Deconvolution" 
Version: 1.0 

Edited based on the references [1][2].

This package contains the implementations of the algorithms described in the paper, including PSF (blur kernel) estimation from a single image, Zernike-based physically constrained aberration correction algorithm using the estimated PSF and deconvolution reconstruction algorithms, which are described in detail in the paper: "Computational Adaptive Optics for Fluorescence Microscopy via Sparse Blind Deconvolution".

Please cite our paper if using the code to generate data (e.g., images, tables of processing times, etc.) in uncommercial usage.


## How to use

The code is tested in MATLAB 2024b (64bit) under the Windows 11 64bit version with an Intel i9-13900KF CPU, NVIDIA GeForce RTX 4060 GPU and 64GB RAM.

1. Unpack the package

2. Include the subdirectory in your MATLAB path

3. Run the .m files with the prefix of "Main" to process the example samples.



## Main modules description

1. k_cal.m:

   * Input a blurred image (Blurred) and the kernel size (kx, ky), estimates the initial PSF from a single blurred image under a sparse prior constraint.

2. Main.m:

   * Performs blind PSF estimation combined with Zernike-based physically constrained aberration correction for deconvolution reconstruction.

## Citation

If you use this code, please cite this paper as well as the paper referred to in this paper:

[1]. Zhang R, Du H, Zhou N, et al. Computational Adaptive Optics for Fluorescence Microscopy via Sparse Blind Deconvolution[J]. Laser & Photonics Reviews, 2025: 2500032. https://doi.org/10.1002/lpor.202500032

[2]. Levin A, Weiss Y, Durand F, et al. Understanding and evaluating blind deconvolution algorithms[C]//2009 IEEE conference on computer vision and pattern recognition. IEEE, 2009: 1964-1971.&#x20; https://doi.org/10.1109/CVPR.2009.5206815

