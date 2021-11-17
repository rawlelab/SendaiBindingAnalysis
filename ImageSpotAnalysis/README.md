Dual Color Analysis

Original version of script:
By Bob Rawle, Kasson Lab, University of Virginia, 2016
Published online in conjunction with:
Rawle et al., Disentangling Viral Membrane Fusion from Receptor Binding Using Synthetic DNA-Lipid Conjugates, Biophysical Journal (2016) 
http://dx.doi.org/10.1016/j.bpj.2016.05.048

This is a slightly modified version of the original code to analyze Sendai virus binding data. The main differences are different parameters and options used in the analysis that are optimized for the Sendai virus data.

Modifications by:
Bob Rawle and Amy Lam, Williams College, 2021
Published online in conjunction with:
Lam et al, 2021, Single virus assay reveals membrane determinants and mechanistic features of Sendai virus binding

To begin:
To start the program, run the function Start_Dual_Color_Analysis.
Before starting the program, the options should be specified in Setup_Options.

Basic description:
This program will analyze a stack of micrograph images of fluorescently labeled particles that have been imaged in two colors. In the case of the Biophysical Journal article, this referred to fluorescently labeled viral particles that also contained fluorescently tagged DNA-lipids. The goal of the program is to find the particles in one color, and then quantify the fluorescence intensity in both colors within the region of interest around each particle. The micrographs should be formatted as a stack of TIF images, 16-bit monochromatic. The default ordering should be the first color image (used to find the particles) followed by the second color image.

The general algorithm is as follows: 
1) The program loads each image in the stack which will be used to locate the fluorescently labeled particles. This image is displayed as Figure 1. The second color image is displayed in Figure 2.
2) A global threshold is applied to generate a logical image, displayed as Figure 3.
3) The logical image is used to identify fluorescently labeled particles.  Each particle is designated as either "good" or "bad", according to a variety of tests – such as whether the particle happens to be too close to the edge of the field of view, too close to another particle, etc. "Good" particles are drawn with a green box in Figure 1, and "bad" particles are drawn with a red box.
4) For each particle, the intensity within the region of interest is calculated for both the first and second color images.
6) The data is saved as a .mat file.

For additional details, see the notes throughout the program scripts.
