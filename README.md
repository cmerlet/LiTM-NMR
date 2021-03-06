# LiTM-NMR

* Program calculating model NMR spectra for Li ions in LiMnxTi(2-x)O4 according to Mn-O-Li pathways obtained through DFT and mean field evaluations of scaling factors. This program can be used for x = 0.5 to x = 2.0.

* To compile the program, you can use any c compiler, e.g.:

      gcc -o LiTM-NMR.x LiTM-NMR.c -lm

* To execute the program, you just have to do:

      ./program < inputfile

* The input file contains the following information:

1. Name of the xyz file with the initial positions of a P4332 structure of LiTi1.5Mn0.54 

2. Length of the box in all dimensions 

3. Number of periodic images in all dimensions 

4. Seed for the random operations, and 0 or 1 for unconstrained/constrained

      0 for unconstrained
      
      Ea in eV for constrained (constrained = no diamagnetic Li, no negative shifts, no Li with a Mn3+ in oct)

5. Fraction of Mn in the system

6. Chemical shifts you want to use

      e.g. HYB20, HYB35, AVDFT, AVmodif, AVdynam
      
      Note that AVdynam only makes sense for x = 1.0 to x = 2.0.

7. Gaussian widths for the model NMR spectrum

      Note that for x > 1.0 and x <= 2.0, the two phases are Fd3m

8. Window and step to plot the NMR spectrum

9. Fraction of Ti-rich region (P4332 region) and fraction of Mn in this region

10. Fraction of Mn tet sites and fraction of Ti in tet site for the first region

      Note that for Mn > 1.0, it is believed that only Ti4+ goes into tet sites, so the 1st value should probably be 0
      
      Note that the fraction of Li in oct sites is the sum of these two fractions

11. Fraction of Mn in tet sites and fraction of Ti in tet site for the second region

      Note that for Mn > 1.0, it is believed that only Ti4+ goes into tet sites, so the 1st value should probably be 0
      
      Note that the fraction of Li in oct sites is the sum of these two fractions

12. Number of steps of additional swaps

      Nequil1 is for Lioct/MnOct/Tioct swaps before the main swaps.
      
      This allows one to switch from P4332 to Fd3m in the Ti-poor region.
      
      Note that for Mn > 1.0, both regions are Fd3m by construction.
      
      Nequil2 is for Lioct/MnOct/Tioct swaps after the main swaps.
      
      Nequil3 is for Lioct/Tioct swaps after the main swaps. 

13. Fraction of Mn3+ replaced by Mn2+/Mn4+ in oct sites

      For now, this is done independently of the Ti-rich / Ti-poor regions.
