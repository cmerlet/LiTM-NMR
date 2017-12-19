# LiTM-NMR

* Program calculating model NMR spectra for Li ions in LiMn(x)Ti(2-x)O4 according to Mn-O-Li pathways obtained through DFT and mean field evaluations of scaling factors. This program can be used for x = 0.5 to x = 2.0.

* To compile the program, you can use any c compiler, e.g.:

`gcc -o LiTM-NMR.x LiTM-NMR.c -lm`

* To execute the program, you have to launch it with an input file giving the following 
information (./program < inputfile):

1. name of the xyz file with the initial positions of a P4332 structure of LiTi1.5Mn0.54 (> name)

2. length of the box in all dimensions (> Lx, Ly, Lz)

3. number of periodic images in all dimensions (> Nx, Ny, Nz)

4. seed for the random operations and unconstrained/constrained (> seed, Ea_const)

0 for unconstrained, Ea in eV for constrained (constrained = no diamagnetic Li, no negative shifts, no Li with a Mn3+ in oct)

5. fraction of Mn in the system (> xMn)

6. chemical shifts you want to use (> shifttype) 
e.g. HYB20, HYB35, AVDFT, AVmodif, AVdynam 
--- Note that AVdynam only makes sense for x=1.0 to x=2.0. ---

7. gaussian widths for the model NMR spectrum (> gwidth for P4332, gwidth2 for Fd3m)
Note that for x>1.0 and x<=2.0, the two phases are Fd3m

8. window and step to plot the NMR spectrum (> smin, smax, dshift2)

9. fraction of Ti-rich region (P4332 region) and fraction of Mn in this region (> pfrac, pxMn)

10. fraction of Mn tet sites and fraction of Ti in tet site for the first region (> invfrac1, invfracTi1)
--- Note that for Mn>1.0, it is believed that only Ti4+ goes into tet sites, so 1st should probably be 0 ---
--- Note that the fraction of Li in oct sites is the sum of these two fractions ---

11. fraction of Mn in tet sites and fraction of Ti in tet site for the second region (> invfrac2, invfracTi2)
--- Note that for Mn>1.0, it is believed that only Ti4+ goes into tet sites, so 1st should probably be 0 ---
--- Note that the fraction of Li in oct sites is the sum of these two fractions ---

12. number of steps of additional swaps (> Nequil1, Nequil2, Nequil3) 
Nequil1 is for Lioct/MnOct/Tioct swaps before the main swaps
--- This allows one to switch from P4332 to Fd3m in the Ti-poor region ---
--- Note that for Mn>1.0, both regions are Fd3m by construction ---
Nequil2 is for Lioct/MnOct/Tioct swaps after the main swaps
Nequil3 is for Lioct/Tioct swaps after the main swaps 

13. fraction of Mn3+ replaced by Mn2+/Mn4+ in oct sites (> oxredfrac)
For now, this is done independently of the Ti-rich / Ti-poor regions.
