/* Program calculating model NMR spectra for Li ions in LiMn(x)Ti(2-x)O4 according to Mn-O-Li pathways obtained through DFT 
and mean field evaluations of scaling factors */
/* This program can be used for x = 0.5 to x = 2.0 */ 

/* To execute the program, you have to launch it with an input file giving the following 
information (./program < inputfile):
1st line: name of the xyz file with the initial positions of a P4332 structure of LiTi1.5Mn0.54 (> name)
2nd line: length of the box in all dimensions (> Lx, Ly, Lz)
3rd line: number of periodic images in all dimensions (> Nx, Ny, Nz)
4th line: seed for the random operations and unconstrained/constrained (> seed, Ea_const)
0 for unconstrained, Ea in eV for constrained (constrained = no diamagnetic Li, no negative shifts, no Li with a Mn3+ in oct)
5th line: fraction of Mn in the system (> xMn)
6th line: chemical shifts you want to use (> shifttype) 
e.g. HYB20, HYB35, AVDFT, AVmodif, AVdynam 
--- Note that AVdynam only makes sense for x=1.0 to x=2.0. ---
7th line: gaussian widths for the model NMR spectrum (> gwidth for P4332, gwidth2 for Fd3m)
Note that for x>1.0 and x<=2.0, the two phases are Fd3m
8th line: window and step to plot the NMR spectrum (> smin, smax, dshift2)
9th line: fraction of Ti-rich region (P4332 region) and fraction of Mn in this region (> pfrac, pxMn)
10th line: fraction of Mn tet sites and fraction of Ti in tet site for the first region (> invfrac1, invfracTi1)
--- Note that for Mn>1.0, it is believed that only Ti4+ goes into tet sites, so 1st should probably be 0 ---
--- Note that the fraction of Li in oct sites is the sum of these two fractions ---
11th line: fraction of Mn in tet sites and fraction of Ti in tet site for the second region (> invfrac2, invfracTi2)
--- Note that for Mn>1.0, it is believed that only Ti4+ goes into tet sites, so 1st should probably be 0 ---
--- Note that the fraction of Li in oct sites is the sum of these two fractions ---
12th line: number of steps of additional swaps (> Nequil1, Nequil2, Nequil3) 
Nequil1 is for Lioct/MnOct/Tioct swaps before the main swaps
--- This allows one to switch from P4332 to Fd3m in the Ti-poor region ---
--- Note that for Mn>1.0, both regions are Fd3m by construction ---
Nequil2 is for Lioct/MnOct/Tioct swaps after the main swaps
Nequil3 is for Lioct/Tioct swaps after the main swaps 
13th line: fraction of Mn3+ replaced by Mn2+/Mn4+ in oct sites (> oxredfrac)
For now, this is done independently of the Ti-rich / Ti-poor regions. */

/* Comment: in bboxtype matrix, -2 is Mn4+, -1 is Mn3+, 0 is Mn2+, 1 is Li+, 2 is Ti4+ and 3 is O2-. */

/* Links to other necessary C files */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
/* L is needed to read input data in external files */
#define L 1000
#define lim 200000
/* Boltzmann constant in eV.K-1 */
#define kB 8.6173324e-5
/* PI number */
#define PI 3.1415926536
/* Constant temperature for now (in K) */
#define Temp 320
/* Rcut ONLY VALID FOR LiTi1.5Mn0.5O4 P4332 STRUCTURE */
#define Rcut 4.0
#define nbins 5000
/* Spinel structure, 3 types of sites, tet/oct/oxygen */
#define Ntypes 3

/* Definition of global variables */ 
/* The default initial values of these variables are zero. */
/* The star(s) define(s) the dimensionality of the numbers/vectors/matrices, these 
entities are then dynamically allocated */
char name[50],shifttype[50],strtype[11][4];
int Ntot,Nx,Ny,Nz,nbins2,Nequil1,Nequil2,Nequil3,seed,startoct,endoct,starttet,endtet,reg1_NMn3,reg1_NMn4,reg2_NMn3,reg2_NMn4;
int *sboxtype,*Ntottype,*counttype,*sitetypeoct;
int **bboxtype,**neighbours_aroundtet,**neighbours_aroundoct;
double Litet_Mn2oct,Lioct_Mn2tet,Litet_Mn3oct,Litet_Mn4oct,Lioct_Mn2oct,Lioct_Mn3oct,Lioct_Mn4oct;
double Lx,Ly,Lz,invfrac1,invfrac2,invfracTi1,invfracTi2,xMn,oxredfrac;
double gwidth,gwidth2,smin,smax,dshift2,Etrap,Mn3frac,pfrac,pxMn,fxMn,ffrac,fXMn,Ea_const,Pacc;
double tMn3frac,pMn3frac,fMn3frac,tTifrac,pTifrac,fTifrac;
double *posinitx,*posinity,*posinitz;
double **posx,**posy,**posz;

/* Specific functions written for this program */
void read();
void find_neighbours();
int count_LioctMn3();
int count_Lidia();
int count_Lineg();
void fromMnTitoMn();
void replaceTiMn3();
void replaceMn4Ti();
void replaceMn3Mn24();
void swapLioctTioct(int Nwanted);
void shuffle_oct(int Nwanted);
void swapLitetMnoct(double invfrac);
void swapLitetTioct(double invfracTi);
void check_stoichio();
void chemical_shift();

/* Functions for the dynamical allocation (these functions are at the end of the program) */
int **imatrix(int nl, int nc);
double *dvector(int n);
double **dmatrix(int nl, int nc);
int *ivector(int n);
double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE);
double *****fiveddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE, int W_SIZE, int V_SIZE);

/* Main function */
int main(void)
{
 clock_t begin,end;
 FILE *outxyz;
 int i,j,NTif,NTip,NTitot,stoptest;
 double time_spent;

 printf("Hello!\n");

 stoptest=0;
 
 /* Starting the clock to see the execution time of the program */
 begin=clock();

 /* Function to read the input data */
 read();
 
 /**********************************************************************/
 /* Part related to x<=1.0, i.e. LiTi1.5Mn0.5O4 to LiTi1.0Mn1.0O4 */

 if(xMn<=1.0)
   {
    /* If there is an Fd3m fraction, check if the user requests can be achieved */
    /* Both regions, Fd3m and P4332 are expected to have xMn between 0.5 and 1.0 */
    fMn3frac=0.0;
    if(ffrac>0)
      {
       NTitot=counttype[2]*(1.5-Mn3frac)/2.0;
       NTip=counttype[2]*pfrac*(1.5-pMn3frac)/2.0;
       NTif=NTitot-NTip;
       fMn3frac=1.5-2.0*NTif/(counttype[2]*ffrac);
      }

    if(fMn3frac>0.5) {printf("\nWe stop the calculation here because we cannot achieve this stoichiometry."); stoptest=1;}

    if(fMn3frac<=0.5)
      {
       if(pfrac>0) printf("\nThe stoichiometry of the P4332 phase will be:LiTi%1.2lfMn%1.2lfO4\n",1.5-pMn3frac,0.5+pMn3frac);
       if(ffrac>0) printf("The stoichiometry of the Fd3m phase will be:LiTi%1.2lfMn%1.2lfO4\n",1.5-fMn3frac,0.5+fMn3frac);

       /* Function to find neighbours */
       find_neighbours();
    
       /* Deal with the P4332 part first */
       if(pfrac>0)
         {
          starttet=1;  endtet=counttype[1]*pfrac;
          startoct=1;  endoct=counttype[2]*pfrac;
          tMn3frac=pMn3frac;
          /* Function to exchange Litet and Mnoct */
          if(invfrac1>0.0) swapLitetMnoct(invfrac1);
          /* Function to exchange Litet and Tioct */
          if(invfracTi1>0.0) swapLitetTioct(invfracTi1);
          /* Function to transform some of the Ti4+/Mn2+ into Mn3+ */
          if(Mn3frac!=0.0) {replaceTiMn3();}
          /* Additional swaps */
          if(Nequil2>0) shuffle_oct(Nequil2);
          if(Nequil3>0) swapLioctTioct(Nequil3);
         }

       /* Deal with the Fd3m part then */
       if(ffrac>0)
         {
          starttet=endtet+1;  endtet=counttype[1];
          startoct=endoct+1;  endoct=counttype[2];
          tMn3frac=fMn3frac;
          /* Function to shuffle Lioct/Mnoct/Tioct before the rest (this allows to start from Fd3m and not P4332) */
          if(Nequil1>0) shuffle_oct(Nequil1);
          /* Function to exchange Litet and Mnoct */
          if(invfrac2>0.0) swapLitetMnoct(invfrac2);
          /* Function to exchange Litet and Tioct */
          if(invfracTi2>0.0) swapLitetTioct(invfracTi2);
          /* Function to transform some of the Ti4+/Mn2+ into Mn3+ */
          if(Mn3frac!=0.0) {replaceTiMn3();}
          /* Additional swaps */
          if(Nequil2>0) shuffle_oct(Nequil2);
          if(Nequil3>0) swapLioctTioct(Nequil3);
         }
      }
   }

 /**********************************************************************/
 
 /**********************************************************************/
 /* Part related to x>1.0, i.e. LiTi1.0Mn1.0O4 to LiMn2.0O4 */

 if(xMn>1.0)
   {
    /* Transform LiTi1.5Mn0.5O4 into LiMn2O4 */
    fromMnTitoMn();

    /* Function to find neighbours */
    find_neighbours();

    /* If there is an Fd3m fraction, check if the user requests can be achieved */
    pTifrac=(2.0-pxMn);
    fTifrac=0.0;
    if(ffrac>0)
      {
       NTitot=counttype[2]*(2.0-xMn)/2.0;
       NTip=counttype[2]*pfrac*(2.0-pxMn)/2.0;
       NTif=NTitot-NTip;
       fTifrac=2.0*NTif/(counttype[2]*ffrac);
       fxMn=2.0-fTifrac;
      }
    
    if(fTifrac>=1.0) {printf("\nWe stop the calculation here because we cannot achieve this stoichiometry."); stoptest=1;}

    if((fTifrac<1.0)&&((2.0-xMn)>0))
      {
       if(pfrac>0) printf("\nThe stoichiometry of the Ti-rich phase will be:LiTi%1.2lfMn%1.2lfO4\n",pTifrac,pxMn);
       if(ffrac>0) printf("The stoichiometry of the Ti-poor phase will be:LiTi%1.2lfMn%1.2lfO4\n",fTifrac,fxMn);

       /* Deal with the Ti-rich part first */
       if(pfrac>0)
         {
          starttet=1;  endtet=counttype[1]*pfrac;
          startoct=1;  endoct=counttype[2]*pfrac;
          tTifrac=pTifrac;
          replaceMn4Ti();
          /* Function to exchange Litet and Mnoct */
	  /* Be careful that this is probably unphysical */
          if(invfrac1>0.0) swapLitetMnoct(invfrac1);
          /* Function to exchange Litet and Tioct */
          if(invfracTi1>0.0) swapLitetTioct(invfracTi1);
          /* Additional swaps */
          if(Nequil2>0) shuffle_oct(Nequil2);
          if(Nequil3>0) swapLioctTioct(Nequil3);
         }

       /* Deal with the Ti-rich part then */
       if(ffrac>0)
         {
          starttet=endtet+1;  endtet=counttype[1];
          startoct=endoct+1;  endoct=counttype[2];
          tTifrac=fTifrac;
          replaceMn4Ti();
          /* Function to exchange Litet and Mnoct */
	  /* Be careful that this is probably unphysical */
          if(invfrac2>0.0) swapLitetMnoct(invfrac2);
          /* Function to exchange Litet and Tioct */
          if(invfracTi2>0.0) swapLitetTioct(invfracTi2);
          /* Additional swaps */
          if(Nequil2>0) shuffle_oct(Nequil2);
          if(Nequil3>0) swapLioctTioct(Nequil3);
         }

      }
   }
 
 /**********************************************************************/

 /* This can be done for all x. */

 /* Check if the program was stopped for stoichiometry issues */
 if(stoptest==0)
   {
    /* Replace some of the Mn3+ by Mn4+ and Mn2+ to see the effect */
    if(oxredfrac>0) replaceMn3Mn24();

    /* Function to check the stoichiometry */ 
    check_stoichio();

    /* Function to calculate the chemical shift of Li and plot the corresponding histogram / NMR spectrum */
    chemical_shift();

    /* Write down the coordinates of the final matrix */
    outxyz=fopen("Positions.xyz","w");
    fprintf(outxyz,"%d\n",Nx*Ny*Nz*(Ntot));
    fprintf(outxyz,"Created with LiTM_NMR.c\n");
    for(i=1;i<=Ntypes;i++)
       {
        for(j=1;j<=counttype[i];j++)
	   {
	    if(i>0) fprintf(outxyz,"%s	%lf	%lf	%lf\n",strtype[bboxtype[i][j]],posx[i][j],posy[i][j],posz[i][j]);
	    if(i==0)  fprintf(outxyz,"Mn2	%lf	%lf	%lf\n",posx[i][j],posy[i][j],posz[i][j]);
	    if(i==-1)  fprintf(outxyz,"Mn3	%lf	%lf	%lf\n",posx[i][j],posy[i][j],posz[i][j]);
	    if(i==-2)  fprintf(outxyz,"Mn4	%lf	%lf	%lf\n",posx[i][j],posy[i][j],posz[i][j]);
	   }
       }
    fclose(outxyz);
   }

 /* Ending the clock and showing the execution time in the terminal */
 end=clock();
 time_spent=(end-begin)/CLOCKS_PER_SEC;
 printf("\nExecution time: %lf.\n",time_spent);

 return 0;
}


/* Read data and allocate matrices and vectors */
void read()
{
 FILE *in;
 char ligne[L],rtype[4];
 int i,j,k,n,Nmax,NLi,NMn,NO,NTi;
 double newx,newy,newz,dx,dy,dz,dist;

 /* Reading of the input file */
 scanf("%s",name);					/* Name of the xyz file with the initial positions, P4332 structure of LiTi1.5Mn0.5O4 */
 scanf("%lf %lf %lf",&Lx,&Ly,&Lz);			/* Length of the box in the three dimensions */
 scanf("%d %d %d",&Nx,&Ny,&Nz);				/* Number of periodic images in three dimensions */
 scanf("%d %lf",&seed,&Ea_const);			/* Seed for the random operations and Ea in eV for acceptance probability of constraints */
 Pacc=exp(-1.0*Ea_const/(kB*Temp));
 scanf("%lf",&xMn);					/* Total fraction of Mn in the system (xMn2 + xMn3 + xMn4) */ 
 if(xMn<=1.0) Mn3frac=xMn-0.5;				/* Fraction for the replacement of Mn2+ and Ti4+ by Mn3+ in the system */
 /* if xMn>1.0, Mn3frac not used so can stay equal to 0 */
 /* Atom types already decided: 0 = Mn, 1 = Li, 2 = Ti, 3 = O */
 strcpy(strtype[0],"Mn");
 strcpy(strtype[1],"Li");
 strcpy(strtype[2],"Ti");
 strcpy(strtype[3],"O");
 scanf("%s",shifttype);					/* Chemical shifts from HYB20, HYB35, AVDFT, AVmodif or AVdynam */
 scanf("%lf %lf",&gwidth,&gwidth2);			/* Gaussian width for the NMR model, gwidth for P4332, gwidth2 for Fd3m */
 scanf("%lf %lf %lf",&smin,&smax,&dshift2);		/* Window and step to plot the NMR spectrum */
 nbins2=(smax-smin)/dshift2;
 scanf("%lf %lf",&pfrac,&pxMn);				/* Fraction of Ti-rich region and xMn in this region */
 if(xMn<=1.0) pMn3frac=pxMn-0.5;
 /* if xMn>1.0, pMn3frac not used so can stay equal to 0 */
 ffrac=1.0-pfrac;
 scanf("%lf %lf",&invfrac1,&invfracTi1);		/* Fraction of Mn in tet sites and fraction of Ti in tet site for the first region */
 scanf("%lf %lf",&invfrac2,&invfracTi2);		/* Fraction of Mn in tet sites and fraction of Ti in tet site for the second region */
 scanf("%d %d %d",&Nequil1,&Nequil2,&Nequil3);		/* Number of steps of reorganisation in oct sites */
 scanf("%lf",&oxredfrac);				/* Fraction of Mn3+ replaced by Mn4+ and Mn2+ in oct sites */

 /* Reading the initial positions, xyz file */
 in=fopen(name,"r");
 fgets(ligne,L,in);
 sscanf(ligne,"%d",&Ntot);				/* Total number of ions/atoms */
 fgets(ligne,L,in);					/* Nothing interesting in this line */
 sboxtype=ivector(Ntot+1);
 Ntottype=ivector(Ntypes+1);
 counttype=ivector(Ntypes+1);
 posinitx=dvector(Ntot+1);
 posinity=dvector(Ntot+1);
 posinitz=dvector(Ntot+1);
 /* Reading atom types and x, y, z positions for all atoms */
 for(i=1;i<=Ntot;i++)
    {
     fgets(ligne,L,in);
     sscanf(ligne,"%s %lf %lf %lf",rtype,&posinitx[i],&posinity[i],&posinitz[i]);
     for(j=0;j<=Ntypes;j++)
	{
	 /* Store the atom type in the sboxtype matrix and count the total number of atoms for each type */
	 if(strcmp(rtype,strtype[j])==0) {sboxtype[i]=j; Ntottype[j]++;}
	}     
    } 
 fclose(in);

 /* Print the total number of atoms for each type and find the max number of atoms, here it is oxygen */
 printf("\nBefore periodic replications:\n");
 Nmax=Ntottype[1];
 for(i=0;i<=Ntypes;i++)
    {
     printf("%s	%d\n",strtype[i],Ntottype[i]);
     if(Ntottype[i]>Nmax) Nmax=Ntottype[i];
    }
 posx=dmatrix(Ntypes+1,Nmax*Nx*Ny*Nz+1);
 posy=dmatrix(Ntypes+1,Nmax*Nx*Ny*Nz+1);
 posz=dmatrix(Ntypes+1,Nmax*Nx*Ny*Nz+1);
 bboxtype=imatrix(Ntypes+1,Nmax*Nx*Ny*Nz+1);
 sitetypeoct=ivector(Nmax*Nx*Ny*Nz+1);

 /* Replicate the box in three dimensions */ 
 NLi=0; NMn=0; NTi=0; NO=0; 
 for(i=1;i<=Nx;i++)
    {
     for(j=1;j<=Ny;j++)
	{
	 for(k=1;k<=Nz;k++)
	    {
	     for(n=1;n<=Ntot;n++)
		{
		 /* Li are all in tetrahedral sites at this point */
		 if(sboxtype[n]==1) 
		   {
		    NLi++;
		    counttype[1]++;
		    newx=posinitx[n]+(i-1)*Lx;
		    newy=posinity[n]+(j-1)*Ly;
		    newz=posinitz[n]+(k-1)*Lz;
		    posx[sboxtype[n]][counttype[sboxtype[n]]]=newx;
		    posy[sboxtype[n]][counttype[sboxtype[n]]]=newy;
		    posz[sboxtype[n]][counttype[sboxtype[n]]]=newz;
		    bboxtype[sboxtype[n]][counttype[sboxtype[n]]]=sboxtype[n];
		   }
		 /* Ti and Mn are all in octahedral sites at this point */
		 if((sboxtype[n]==2)||(sboxtype[n]==0)) 
		   {
		    if(sboxtype[n]==0) NMn++;
		    if(sboxtype[n]==2) NTi++;
		    counttype[2]++;
		    newx=posinitx[n]+(i-1)*Lx;
		    newy=posinity[n]+(j-1)*Ly;
		    newz=posinitz[n]+(k-1)*Lz;
		    posx[2][counttype[2]]=newx;
		    posy[2][counttype[2]]=newy;
		    posz[2][counttype[2]]=newz;
		    bboxtype[2][counttype[2]]=sboxtype[n];
		    if((n<=20)&&(n>=9)) sitetypeoct[counttype[2]]=12;  
		    if((n<=24)&&(n>=21)) sitetypeoct[counttype[2]]=4;  
		   }
		 /* O are only oxygen sites */
		 if(sboxtype[n]==3) 
		   {
		    NO++;
		    counttype[3]++;
		    newx=posinitx[n]+(i-1)*Lx;
		    newy=posinity[n]+(j-1)*Ly;
		    newz=posinitz[n]+(k-1)*Lz;
		    posx[sboxtype[n]][counttype[sboxtype[n]]]=newx;
		    posy[sboxtype[n]][counttype[sboxtype[n]]]=newy;
		    posz[sboxtype[n]][counttype[sboxtype[n]]]=newz;
		    bboxtype[sboxtype[n]][counttype[sboxtype[n]]]=sboxtype[n];
		   }
		}
	    }
	}
    }
 printf("\nAfter periodic replications:\n");
 printf("Total number of tetrahedral sites: %d\n",counttype[1]);
 printf("Total number of Li: %d\n",NLi);
 printf("Total number of octahedral sites: %d\n",counttype[2]);
 printf("Total number of Ti: %d\n",NTi);
 printf("Total number of Mn: %d\n",NMn);
 printf("Total number of oxygen: %d\n",counttype[3]);
 /* Need to change Lx, Ly, Lz because the box as become much bigger */
 Lx=Lx*Nx*1.0;
 Ly=Ly*Ny*1.0;
 Lz=Lz*Nz*1.0;
 
 printf("\nThe acceptance probability of the Monte Carlo moves is %e.\n",Pacc);

 printf("\nReading OK!\n");

}


/* Function to check the actual stoichiometry of the system and calculate fractional occupancies */
void check_stoichio()
{
 FILE *outpdb;
 int i,j,k,NLitet,NLioct,NMn2tet,NMn2oct,NTitet,NTioct,NMn3oct,NMn4oct;
 int ndifftet,ndiffoct,counttet,countoct;
 int *occLitet,*occLioct,*occMn2tet,*occMn2oct,*occTitet,*occTioct,*occMn3oct,*occMn4oct;
 double xLitet,xLioct,xMn2tet,xMn2oct,xMn3oct,xMn4oct,xTitet,xTioct,tempx,tempy,tempz;
 double site8c[6],site12d[6],site4b[6],site8c_pfrac[6],site8c_ffrac[6],site12d_pfrac[6],site12d_ffrac[6],site4b_pfrac[6],site4b_ffrac[6];
 double site16d[6];

 for(i=1;i<=5;i++) 
    {
     site8c[i]=0.0; site12d[i]=0.0; site4b[i]=0.0;
     site8c_pfrac[i]=0.0; site12d_pfrac[i]=0.0; site4b_pfrac[i]=0.0;
     site8c_ffrac[i]=0.0; site12d_ffrac[i]=0.0; site4b_ffrac[i]=0.0;
    }

 NLitet=0; NLioct=0; NTitet=0; NTioct=0; NMn2tet=0; NMn2oct=0; NMn3oct=0; NMn4oct=0; 
 ndifftet=counttype[1]/(Nx*Ny*Nz); ndiffoct=counttype[2]/(Nx*Ny*Nz); 
 occLitet=ivector(ndifftet+1); occMn2tet=ivector(ndifftet+1); occTitet=ivector(ndifftet+1); occTioct=ivector(ndiffoct+1);
 occLioct=ivector(ndiffoct+1); occMn2oct=ivector(ndiffoct+1); occMn3oct=ivector(ndiffoct+1); occMn4oct=ivector(ndiffoct+1);
 for(i=1;i<=ndifftet;i++) {occLitet[i]=0; occMn2tet[i]=0; occTitet[i]=0;}
 for(i=1;i<=ndiffoct;i++) {occLioct[i]=0; occMn2oct[i]=0; occMn3oct[i]=0; occMn4oct[i]=0; occTioct[i]=0;}
 counttet=1; countoct=1;
 for(i=1;i<=counttype[1];i++)
    {
     if(bboxtype[1][i]==0) 
       {
        NMn2tet++; occMn2tet[counttet]++; counttet++;
	if(i<=counttype[1]*pfrac) site8c_pfrac[2]++;
	if(i>counttype[1]*pfrac) site8c_ffrac[2]++;
       }
     if(bboxtype[1][i]==1) 
       {
	NLitet++; occLitet[counttet]++; counttet++;
	if(i<=counttype[1]*pfrac) site8c_pfrac[1]++;
	if(i>counttype[1]*pfrac) site8c_ffrac[1]++;
       }
     if(bboxtype[1][i]==2) 
       {
	NTitet++; occTitet[counttet]++; counttet++;
	if(i<=counttype[1]*pfrac) site8c_pfrac[3]++;
	if(i>counttype[1]*pfrac) site8c_ffrac[3]++;
       }
     if(counttet>ndifftet) counttet=1;
    }
 for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==-2) 
       {
	NMn4oct++; occMn4oct[countoct]++; countoct++;
	if(sitetypeoct[i]==12) 
	  {
	   site12d[4]++;
	   if(i<=counttype[2]*pfrac) site12d_pfrac[4]++;
	   if(i>counttype[2]*pfrac) site12d_ffrac[4]++;
	  }  
	if(sitetypeoct[i]==4) 
	  {
	   site4b[4]++;  
	   if(i<=counttype[2]*pfrac) site4b_pfrac[4]++;
	   if(i>counttype[2]*pfrac) site4b_ffrac[4]++;
	  }
	if(i<=counttype[2]*pfrac) reg1_NMn4++;
	if(i>counttype[2]*pfrac) reg2_NMn4++;
       }
     if(bboxtype[2][i]==-1) 
       {
	NMn3oct++; occMn3oct[countoct]++; countoct++;
	if(sitetypeoct[i]==12) 
	  {
	   site12d[3]++;
	   if(i<=counttype[2]*pfrac) site12d_pfrac[3]++;
	   if(i>counttype[2]*pfrac) site12d_ffrac[3]++;
	  }  
	if(sitetypeoct[i]==4) 
	  {
	   site4b[3]++;
	   if(i<=counttype[2]*pfrac) site4b_pfrac[3]++;
	   if(i>counttype[2]*pfrac) site4b_ffrac[3]++;
	  }  
	if(i<=counttype[2]*pfrac) reg1_NMn3++;
	if(i>counttype[2]*pfrac) reg2_NMn3++;
       }
     if(bboxtype[2][i]==0) 
       {
	NMn2oct++; occMn2oct[countoct]++; countoct++;
	if(sitetypeoct[i]==12) 
	  {
	   site12d[2]++;
	   if(i<=counttype[2]*pfrac) site12d_pfrac[2]++;
	   if(i>counttype[2]*pfrac) site12d_ffrac[2]++;
	  }  
	if(sitetypeoct[i]==4) 
	  {
	   site4b[2]++;
	   if(i<=counttype[2]*pfrac) site4b_pfrac[2]++;
	   if(i>counttype[2]*pfrac) site4b_ffrac[2]++;
	  }  
       }
     if(bboxtype[2][i]==1) 
       {
	NLioct++; occLioct[countoct]++; countoct++;
	if(sitetypeoct[i]==12) 
	  {
	   site12d[1]++;
	   if(i<=counttype[2]*pfrac) site12d_pfrac[1]++;
	   if(i>counttype[2]*pfrac) site12d_ffrac[1]++;
	  }  
	if(sitetypeoct[i]==4) 
	  {
	   site4b[1]++;
	   if(i<=counttype[2]*pfrac) site4b_pfrac[1]++;
	   if(i>counttype[2]*pfrac) site4b_ffrac[1]++;
	  }  
       }
     if(bboxtype[2][i]==2) 
       {
	NTioct++; occTioct[countoct]++; countoct++;
	if(sitetypeoct[i]==12) 
	  {
	   site12d[5]++;
	   if(i<=counttype[2]*pfrac) site12d_pfrac[5]++;
	   if(i>counttype[2]*pfrac) site12d_ffrac[5]++;
	  }  
	if(sitetypeoct[i]==4) 
	  {
	   site4b[5]++; 
	   if(i<=counttype[2]*pfrac) site4b_pfrac[5]++;
	   if(i>counttype[2]*pfrac) site4b_ffrac[5]++;
	  } 
       }
     if(countoct>ndiffoct) countoct=1;
    }
 xLitet=NLitet*1.0/(counttype[1]*1.0);
 xMn2tet=NMn2tet*1.0/(counttype[1]*1.0);
 xTitet=NTitet*1.0/(counttype[1]*1.0);
 xLioct=NLioct*1.0/(counttype[2]*1.0);
 xMn2oct=NMn2oct*1.0/(counttype[2]*1.0);
 xMn3oct=NMn3oct*1.0/(counttype[2]*1.0);
 xMn4oct=NMn4oct*1.0/(counttype[2]*1.0);
 xTioct=NTioct*1.0/(counttype[2]*1.0);
 printf("The system is:");
 printf(" Li %lf Ti %lf Mn %lf [ Li %lf Ti %lf Mn2+ %lf Mn3+ %lf Mn4+ %lf ] O4\n",xLitet,xTitet,xMn2tet,xLioct*2.0,xTioct*2.0,xMn2oct*2.0,xMn3oct*2.0,xMn4oct*2.0);
 printf("Note that in tet, only Mn2+ are allowed\n");

 /* Write a PDB file with fractional occupancies */
 /* For now, only for P4332 fraction when xMn<=1.0 */
  if(xMn<=1.0)
   {
    outpdb=fopen("Initial_box_fractional_occupancies.pdb","w");
    fprintf(outpdb,"CRYST1    %lf    %lf    %lf  90.00  90.00  90.00 P 1           1\n",Lx/(Nx*1.0),Ly/(Ny*1.0),Lz/(Nz*1.0));
    j=1;
    for(i=1;i<=ndifftet;i++)
       {
        tempx=posx[1][i]; tempy=posy[1][i]; tempz=posz[1][i]; 	
        fprintf(outpdb,"ATOM	%d	Li	X 1	%lf	%lf	%lf	%lf	0.00 LI\n",j,tempx,tempy,tempz,occLitet[i]*ndifftet*1.0/(counttype[1]*1.0));
        j++;     
        site8c[1]+=occLitet[i]*ndifftet*1.0/(counttype[1]*1.0);
        fprintf(outpdb,"ATOM	%d	Mn2	X 1	%lf	%lf	%lf	%lf	0.00 MN\n",j,tempx,tempy,tempz,occMn2tet[i]*ndifftet*1.0/(counttype[1]*1.0));
        j++;
        site8c[2]+=occMn2tet[i]*ndifftet*1.0/(counttype[1]*1.0);
        fprintf(outpdb,"ATOM	%d	Ti	X 1	%lf	%lf	%lf	%lf	0.00 Ti\n",j,tempx,tempy,tempz,occTitet[i]*ndifftet*1.0/(counttype[1]*1.0));
        j++;
        site8c[3]+=occTitet[i]*ndifftet*1.0/(counttype[1]*1.0);
       }
    for(i=1;i<=3;i++) {site8c[i]/=8.0; site8c_pfrac[i]/=8.0; site8c_ffrac[i]/=8.0;}
    for(i=1;i<=ndiffoct;i++)
       {
        tempx=posx[2][i]; tempy=posy[2][i]; tempz=posz[2][i]; 	
        fprintf(outpdb,"ATOM	%d	Li	X 1	%lf	%lf	%lf	%lf	0.00 LI\n",j,tempx,tempy,tempz,occLioct[i]*ndiffoct*1.0/(counttype[2]*1.0));
        j++;
        fprintf(outpdb,"ATOM	%d	Mn2	X 1	%lf	%lf	%lf	%lf	0.00 MN\n",j,tempx,tempy,tempz,occMn2oct[i]*ndiffoct*1.0/(counttype[2]*1.0));
        j++;
        fprintf(outpdb,"ATOM	%d	Mn3	X 1	%lf	%lf	%lf	%lf	0.00 MN\n",j,tempx,tempy,tempz,occMn3oct[i]*ndiffoct*1.0/(counttype[2]*1.0));
        j++;
        fprintf(outpdb,"ATOM	%d	Ti	X 1	%lf	%lf	%lf	%lf	0.00 Ti\n",j,tempx,tempy,tempz,occTioct[i]*ndiffoct*1.0/(counttype[2]*1.0));
        j++;
       }
    for(i=1;i<=5;i++) {site12d[i]/=12.0; site4b[i]/=4.0; site12d_pfrac[i]/=12.0; site4b_pfrac[i]/=4.0; site12d_ffrac[i]/=12.0; site4b_ffrac[i]/=4.0;} 
    for(i=1;i<=(counttype[3]/(Nx*Ny*Nz));i++)
       {
        tempx=posx[3][i]; tempy=posy[3][i]; tempz=posz[3][i]; 	
        fprintf(outpdb,"ATOM	%d	O	X 1	%lf	%lf	%lf	%lf	0.00 O\n",j,tempx,tempy,tempz,1.0);
        j++;
       }
    fprintf(outpdb,"END\n");
    fclose(outpdb);
    printf("\nFractional occupancies for the full structure:\n");
    printf("8c  %.3lf Li, %.3lf Mn2+, %.3lf Ti4+:\n",site8c[1],site8c[2],site8c[3]);
    for(i=1;i<=3;i++) {site8c_pfrac[i]*=ndifftet*1.0/(counttype[1]*pfrac*1.0); site8c_ffrac[i]*=ndifftet*1.0/(counttype[1]*ffrac*1.0);}
    for(i=1;i<=5;i++) {site12d[i]*=ndiffoct*1.0/(counttype[2]*1.0); site12d_pfrac[i]*=ndiffoct*1.0/(counttype[2]*pfrac*1.0); 
		       site12d_ffrac[i]*=ndiffoct*1.0/(counttype[2]*ffrac*1.0);} 
    for(i=1;i<=5;i++) {site4b[i]*=ndiffoct*1.0/(counttype[2]*1.0); site4b_pfrac[i]*=ndiffoct*1.0/(counttype[2]*pfrac*1.0); 
		       site4b_ffrac[i]*=ndiffoct*1.0/(counttype[2]*ffrac*1.0);} 
    printf("4b  %.3lf Li, %.3lf Mn2+, %.3lf Mn3+, %.3lf Ti4+:\n",site4b[1],site4b[2],site4b[3],site4b[5]);
    printf("12d  %.3lf Li, %.3lf Mn2+, %.3lf Mn3+, %.3lf Ti4+:\n",site12d[1],site12d[2],site12d[3],site12d[5]);
    if(pfrac>0)
      {
       printf("Fractional occupancies for the P4332 phase:\n");
       printf("8c  %.3lf Li, %.3lf Mn2+, %.3lf Ti4+:\n",site8c_pfrac[1],site8c_pfrac[2],site8c_pfrac[3]);
       printf("4b  %.3lf Li, %.3lf Mn2+, %.3lf Mn3+, %.3lf Ti4+:\n",site4b_pfrac[1],site4b_pfrac[2],site4b_pfrac[3],site4b_pfrac[5]);
       printf("12d  %.3lf Li, %.3lf Mn2+, %.3lf Mn3+, %.3lf Ti4+:\n",site12d_pfrac[1],site12d_pfrac[2],site12d_pfrac[3],site12d_pfrac[5]);
      }
    if(ffrac>0)
      {
       printf("Fractional occupancies for the Fd3m phase:\n");
       for(i=1;i<=5;i++) site16d[i]=(site4b_ffrac[i]+site12d_ffrac[i])/2.0;
       printf("8a  %.3lf Li, %.3lf Mn2+, %.3lf Ti4+:\n",site8c_ffrac[1],site8c_ffrac[2],site8c_ffrac[3]);
       printf("16d  %.3lf Li, %.3lf Mn2+, %.3lf Mn3+, %.3lf Ti4+:\n",site16d[1],site16d[2],site16d[3],site16d[5]);
      }
   }
 

}


/* Function to find neighbours around a tetrahedral or an octahedral site */
void find_neighbours()
{
 int i,j,Nneigh;
 double dx,dy,dz,dist;
 
 neighbours_aroundtet=imatrix(counttype[1]+1,13);
 neighbours_aroundoct=imatrix(counttype[2]+1,13);

 /* Loop over the octahedral sites */
 for(i=1;i<=counttype[2];i++)
    {
     Nneigh=0;
     /* Find the neighbours in oct positions */
     for(j=1;j<=counttype[2];j++)
	{
	 if(j!=i)
	   {
	    dx=posx[2][i]-posx[2][j];
	    dy=posy[2][i]-posy[2][j];
	    dz=posz[2][i]-posz[2][j];
            /* Periodic boundary conditions */
            if(dx>(Lx/2.0)) dx-=Lx;
            if(dx<(-Lx/2.0)) dx+=Lx;
            if(dy>(Ly/2.0)) dy-=Ly;
            if(dy<(-Ly/2.0)) dy+=Ly;
            if(dz>(Lz/2.0)) dz-=Lz;
            if(dz<(-Lz/2.0)) dz+=Lz;
	    dist=dx*dx+dy*dy+dz*dz;
	    dist=sqrt(dist);
	    if(dist<=Rcut) {Nneigh++; neighbours_aroundoct[i][Nneigh]=j;}
           }
        }
     /* Find the neighbours in tet positions */
     for(j=1;j<=counttype[1];j++)
	{
	 dx=posx[2][i]-posx[1][j];
	 dy=posy[2][i]-posy[1][j];
	 dz=posz[2][i]-posz[1][j];
         /* Periodic boundary conditions */
         if(dx>(Lx/2.0)) dx-=Lx;
         if(dx<(-Lx/2.0)) dx+=Lx;
         if(dy>(Ly/2.0)) dy-=Ly;
         if(dy<(-Ly/2.0)) dy+=Ly;
         if(dz>(Lz/2.0)) dz-=Lz;
         if(dz<(-Lz/2.0)) dz+=Lz;
	 dist=dx*dx+dy*dy+dz*dz;
	 dist=sqrt(dist);
	 if(dist<=Rcut) {Nneigh++; neighbours_aroundoct[i][Nneigh]=j;}
        }
     if(Nneigh!=12) printf("Problem with %d, Nneigh (oct) = %d\n",i,Nneigh); 
    }

 /* Loop over the tetrahedral sites */
 for(i=1;i<=counttype[1];i++)
    {
     Nneigh=0;
     /* Find the neighbours in oct positions */
     for(j=1;j<=counttype[2];j++)
        {
         dx=posx[1][i]-posx[2][j];
         dy=posy[1][i]-posy[2][j];
         dz=posz[1][i]-posz[2][j];
         /* Periodic boundary conditions */
         if(dx>(Lx/2.0)) dx-=Lx; 
         if(dx<(-Lx/2.0)) dx+=Lx; 
         if(dy>(Ly/2.0)) dy-=Ly; 
         if(dy<(-Ly/2.0)) dy+=Ly; 
         if(dz>(Lz/2.0)) dz-=Lz; 
         if(dz<(-Lz/2.0)) dz+=Lz; 
         dist=dx*dx+dy*dy+dz*dz;
         dist=sqrt(dist);
         if(dist<=Rcut) {Nneigh++; neighbours_aroundtet[i][Nneigh]=j;}
        }
     if(Nneigh!=12) printf("Problem with %d, Nneigh (tet) = %d\n",i,Nneigh); 
    }

}


/* Count the number of Lioct having a Mn3+ in its neighbours */
int count_LioctMn3()
{
 int i,n,NLioctMn3,NMn3;

 /* NMn3 here is the number of Mn3+ in the neighbours */

 NLioctMn3=0;
 /* Check for Li in oct */
 for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==1)
       {
        NMn3=0;
        /* First 6 neighbours are in oct sites */
        for(n=1;n<=6;n++) {if((bboxtype[2][neighbours_aroundoct[i][n]])==-1) NMn3++;}
        /* Last 6 neighbours are in tet sites */
        for(n=7;n<=12;n++) {if((bboxtype[1][neighbours_aroundoct[i][n]])==-1) NMn3++;}
        if(NMn3!=0) NLioctMn3++;
       }
    }

 return NLioctMn3;
}


/* Count the number of diamagnetic lithium in the system */
int count_Lidia()
{
 int i,n,NLidia,NMn;

 /* NMn here is the total number of Mn in the neighbours (Mn2+, Mn3+ and Mn4+) */

 NLidia=0;
 /* Check for Li in tet */
 for(i=1;i<=counttype[1];i++)
    {
     if(bboxtype[1][i]==1)
       {
	NMn=0;
        for(n=1;n<=12;n++) {if((bboxtype[2][neighbours_aroundtet[i][n]])<=0) NMn++;}
        if(NMn==0) NLidia++;
       }
    }
 /* Check for Li in oct */
 for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==1)
       {
	NMn=0;
	/* First 6 neighbours are in oct sites */
        for(n=1;n<=6;n++) {if((bboxtype[2][neighbours_aroundoct[i][n]])<=0) NMn++;}
	/* Last 6 neighbours are in tet sites */
        for(n=7;n<=12;n++) {if((bboxtype[1][neighbours_aroundoct[i][n]])<=0) NMn++;}
        if(NMn==0) NLidia++;
       }
    }
 
 return NLidia;
}


/* Count the number of lithium with a negative shift in the system */
int count_Lineg()
{
 int i,n,NLineg,NMntet,NMn2oct,NMn3oct,NMn4oct,pNMn3oct,fNMn3oct,pNMn4oct,fNMn4oct;
 double avoxstate,pavoxstate,favoxstate,pLioct_Mn3oct,pLioct_Mn4oct,fLioct_Mn3oct,fLioct_Mn4oct;

 /* Shifts (ppm) obtained from HYB20 calculations */
 if(strcmp(shifttype,"HYB20")==0) 
   {Lioct_Mn2tet=18.1; Lioct_Mn2oct=-72.5; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;} 

 /* Shifts (ppm) obtained from HYB35 calculations */
 if(strcmp(shifttype,"HYB35")==0) 
   {Lioct_Mn2tet=15.3; Lioct_Mn2oct=-54.1; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;} 

 /* Shifts (ppm) obtained as average from HYB20 and HYB35 calculations */
 if(strcmp(shifttype,"AVDFT")==0) 
   {Lioct_Mn2tet=16.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;} 
 if(strcmp(shifttype,"AVmodif")==0) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;} 
 /* Shifts (ppm) obtained for dynamic averaging of Mn3+ and Mn4+ in oct */
 /* Here same as AVmodif */
 /* First case, there is actually only one phase */
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac==0)||(ffrac==0))) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4; 
    avoxstate=(4.0*xMn-1.0)/xMn; 
    Lioct_Mn3oct=(4.0-avoxstate)*Lioct_Mn3oct+(avoxstate-3.0)*Lioct_Mn4oct; Lioct_Mn4oct=Lioct_Mn3oct;}
 /* Second case, there are two phases */
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac!=0)&&(ffrac!=0))) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;
    pavoxstate=(4.0*pxMn-1.0)/pxMn; 
    favoxstate=(4.0*fxMn-1.0)/fxMn; 
    pLioct_Mn3oct=(4.0-pavoxstate)*Lioct_Mn3oct+(pavoxstate-3.0)*Lioct_Mn4oct; pLioct_Mn4oct=pLioct_Mn3oct;
    fLioct_Mn3oct=(4.0-favoxstate)*Lioct_Mn3oct+(favoxstate-3.0)*Lioct_Mn4oct; fLioct_Mn4oct=fLioct_Mn3oct;}


 NLineg=0;
 /* Check for Li in oct */
 /* Here I check only for Mn2+ because Mn3+/Mn4+ do not go in tet sites and Lioct-Mn3+oct is not supposed to happen */
 if((strcmp(shifttype,"AVdynam")!=0)||((pfrac==0)||(ffrac==0))) 
   {
    for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==1)
       {
	NMn2oct=0; NMn3oct=0; NMn4oct=0; NMntet=0;
	/* First 6 neighbours are in oct sites */
        for(n=1;n<=6;n++) 
	   {
	    if((bboxtype[2][neighbours_aroundoct[i][n]])==0) NMn2oct++;
	    if((bboxtype[2][neighbours_aroundoct[i][n]])==-1) NMn3oct++;
	    if((bboxtype[2][neighbours_aroundoct[i][n]])==-2) NMn4oct++;
	   }
	/* Last 6 neighbours are in tet sites */
        for(n=7;n<=12;n++) {if((bboxtype[1][neighbours_aroundoct[i][n]])==0) NMntet++;}
        if((NMntet*Lioct_Mn2tet+NMn2oct*Lioct_Mn2oct+NMn3oct*Lioct_Mn3oct+NMn4oct*Lioct_Mn4oct)<0) NLineg++; 
       }
    }
   }
 
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac!=0)&&(ffrac!=0))) 
   {
    for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==1)
       {
	NMn2oct=0; NMn3oct=0; NMn4oct=0; NMntet=0;
	pNMn3oct=0; pNMn4oct=0; fNMn3oct=0; fNMn4oct=0;
	/* First 6 neighbours are in oct sites */
        for(n=1;n<=6;n++) 
	   {
	    if((bboxtype[2][neighbours_aroundoct[i][n]])==0) NMn2oct++;
	    if(((bboxtype[2][neighbours_aroundoct[i][n]])==-1)&&(neighbours_aroundoct[i][n]<=counttype[2]*pfrac)) {NMn3oct++; pNMn3oct++;}
	    if(((bboxtype[2][neighbours_aroundoct[i][n]])==-2)&&(neighbours_aroundoct[i][n]<=counttype[2]*pfrac)) {NMn4oct++; pNMn4oct++;}
	    if(((bboxtype[2][neighbours_aroundoct[i][n]])==-1)&&(neighbours_aroundoct[i][n]>counttype[2]*pfrac)) {NMn3oct++; fNMn3oct++;}
	    if(((bboxtype[2][neighbours_aroundoct[i][n]])==-2)&&(neighbours_aroundoct[i][n]>counttype[2]*pfrac)) {NMn4oct++; fNMn4oct++;}
	   }
	/* Last 6 neighbours are in tet sites */
        for(n=7;n<=12;n++) {if((bboxtype[1][neighbours_aroundoct[i][n]])==0) NMntet++;}
        if((NMntet*Lioct_Mn2tet+NMn2oct*Lioct_Mn2oct+pNMn3oct*pLioct_Mn3oct+pNMn4oct*pLioct_Mn4oct+fNMn3oct*fLioct_Mn3oct+fNMn4oct*fLioct_Mn4oct)<0) NLineg++; 
       }
    }
   }
 
 return NLineg;
}


/* Function to transform all Mn2+/Ti4+ in Mn3+ and then Mn3+/Mn4+ */
/* This is just to go from LiTi1.5Mn0.5O4 to LiMn2O4 */
void fromMnTitoMn()
{
 int i,iMn,Nwanted,Nobs,n,l1,l2;

 printf("Entering fromMnTitoMn\n"); 

 /* First, transform all Ti4+/Mn2+ into Mn3+ */
 for(i=1;i<=counttype[2];i++) bboxtype[2][i]=-1;

 /* Then replace randomly half of the Mn3+ by Mn4+ */
 srand(seed*40+1);
 
 Nwanted=counttype[2]/2;

 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a Mn3 */
     iMn=rand()%counttype[2]+1;
     /* If we are not on a Mn3 site, we move to a Mn3 site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iMn]!=-1)&&(l2<=lim))
          {
     	   iMn=rand()%counttype[2]+1;
           l2++;
          }
     /* Replace this Mn3+ by a Mn4+ */
     if((bboxtype[2][iMn])==-1)
       {
        bboxtype[2][iMn]=-2;
        n++;
       }
    }

 /* Count the actual number of Mn3+ in the octrahedral lattice */
 Nobs=0;
 for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==-2) Nobs++;
    }
 printf("Total number of Mn4+ in the octahedral lattice: %d\n",Nobs);
 printf("This corresponds to a fraction of: %lf\n",Nobs*1.0/(counttype[2]*1.0));

 printf("Leaving fromMnTitoMn\n"); 

}


/* Function to transform some of the Mn4+ into Ti4+ */
void replaceMn4Ti()
{
 int i,l1,l2,n,Nobs,Nwanted,iMn4;
 int NLidia,NLidia_beg,NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg;
 double randnb;

 printf("Entering replace Mn4 Ti\n"); 

 srand(seed+312);

 /* Number of Ti4+ we want in the oct sites */
 Nwanted=(endoct-startoct+1);
 Nwanted*=tTifrac;
 Nwanted/=2;

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the replacement Ti-Mn with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the replacement Ti-Mn with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the replacement Ti-Mn with NLioctMn3 = %d\n",NLioctMn3_beg); 

 /* Replacing Mn4+ by Ti4+ */
 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a Mn4+ */
     iMn4=rand()%(endoct-startoct+1)+startoct;
     /* If we are not on a Mn4+ site, we move to a Mn4+ site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iMn4]!=-2)&&(l2<=lim))
          {
           iMn4=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Replace this Mn4+ by a Ti4+ */
     if((bboxtype[2][iMn4])==-2)
       {
        bboxtype[2][iMn4]=2;
        n++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>=Pacc))
	{
	 bboxtype[2][iMn4]=-2;
	 n--;
	}
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
    }

 /* Count the actual number of Ti4+ in the octrahedral lattice */
 Nobs=0;
 for(i=startoct;i<=endoct;i++)
    {
     if(bboxtype[2][i]==2) Nobs++;
    }
 printf("Total number of Ti4+ in the octahedral lattice: %d\n",Nobs);
 printf("This corresponds to a fraction of: %lf\n",Nobs*1.0/((endoct-startoct+1)*1.0));

 printf("Leaving replace Mn4 Ti\n"); 
}


/* Function to transform some of the Ti4+ and Mn2+ into Mn3+ */
void replaceTiMn3()
{
 int i,l1,l2,n,Nobs,Nwanted,iTi,iMn2;
 int NLidia,NLidia_beg,NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg;
 double randnb;

 printf("Entering replace Ti Mn3\n"); 

 srand(seed*12+14);

 /* Number of Mn3+ we want in the oct sites system */
 Nwanted=(endoct-startoct+1);
 Nwanted*=tMn3frac;
 /* Divide by two to replace Ti4+ and Mn2+ separately */
 Nwanted/=2;

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the replacement Ti-Mn with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the replacement Ti-Mn with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the replacement Ti-Mn with NLioctMn3 = %d\n",NLioctMn3_beg); 

 /* Replacing Ti4+ by Mn3+ */
 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a titanium */
     iTi=rand()%(endoct-startoct+1)+startoct;
     /* If we are not on a Ti site, we move to a Ti site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iTi]!=2)&&(l2<=lim))
          {
           iTi=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Replace this Ti4+ by a Mn3+ */
     if((bboxtype[2][iTi])==2)
       {
        bboxtype[2][iTi]=-1;
        n++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>=Pacc))
	{
	 bboxtype[2][iTi]=2;
	 n--;
	}
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
    }

 /* Replacing Mn2+ by Mn3+ */
 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a Mn2+ */
     iMn2=rand()%(endoct-startoct+1)+startoct;
     /* If we are not on a Mn2+ site, we move to a Mn2+ site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iMn2]!=0)&&(l2<=lim))
          {
           iMn2=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Try to replace this Mn2+ by a Mn3+ */
     if((bboxtype[2][iMn2])==0)
       {
        bboxtype[2][iMn2]=-1;
        n++;
       }
     /* Checking the constraints: maybe this is not necessary because Mn2+/Mn3+ are very similar */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
	{
	 bboxtype[2][iMn2]=0;
	 n--;
	}
    }

 
 /* Count the actual number of Mn3+ in the octrahedral lattice */
 Nobs=0;
 for(i=startoct;i<=endoct;i++)
    {
     if(bboxtype[2][i]==-1) Nobs++;
    }
 printf("Total number of Mn3+ in the octahedral lattice: %d\n",Nobs);
 printf("This corresponds to a fraction of: %lf\n",Nobs*1.0/((endoct-startoct+1)*1.0));

 printf("Leaving replace Ti Mn3\n"); 
}


/* Function to transform some of the Mn3+ into Mn2+ and Mn4+ */
void replaceMn3Mn24()
{
 int i,l1,l2,n,Nobs2,Nobs4,Nwanted,iMn3;
 int NLidia,NLidia_beg,NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg;
 double randnb;

 printf("Entering replace Mn3 by Mn2 Mn4\n"); 

 srand(seed*80+42);

 /* Number of Mn3+ we want to replace in the oct sites system */
 Nwanted=counttype[2]*oxredfrac;
 /* Divide by two to replace Mn4+ and Mn2+ separately */
 Nwanted/=2;

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the replacement Ti-Mn with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the replacement Ti-Mn with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the replacement Ti-Mn with NLioctMn3 = %d\n",NLioctMn3_beg); 

 /* Replacing Mn3+ by Mn4+ */
 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a Mn3+ */
     iMn3=rand()%(counttype[2])+1;
     /* If we are not on a Mn3+ site, we move to a Mn3+ site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iMn3]!=-1)&&(l2<=lim))
          {
           iMn3=rand()%(counttype[2])+1;
           l2++;
          }
     /* Replace this Mn3+ by a Mn4+ */
     if((bboxtype[2][iMn3])==-1)
       {
        bboxtype[2][iMn3]=-2;
        n++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>=Pacc))
	{
	 bboxtype[2][iMn3]=-1;
	 n--;
	}
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
    }

 /* Replacing Mn3+ by Mn2+ */
 n=0;
 l1=1;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* Randomly choose a Mn3+ */
     iMn3=rand()%(counttype[2])+1;
     /* If we are not on a Mn3+ site, we move to a Mn3+ site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iMn3]!=-1)&&(l2<=lim))
          {
           iMn3=rand()%(counttype[2])+1;
           l2++;
          }
     /* Try to replace this Mn3+ by a Mn2+ */
     if((bboxtype[2][iMn3])==-1)
       {
        bboxtype[2][iMn3]=0;
        n++;
       }
     /* Checking the constraints: maybe this is not necessary because Mn2+/Mn3+ are very similar */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
	{
	 bboxtype[2][iMn3]=-1;
	 n--;
	}
    }

 
 /* Count the actual number of Mn3+ in the octrahedral lattice */
 Nobs2=0; Nobs4=0;
 for(i=1;i<=counttype[2];i++)
    {
     if(bboxtype[2][i]==-2) Nobs4++;
     if(bboxtype[2][i]==0) Nobs2++;
    }
 printf("Total number of Mn4+ in the octahedral lattice: %d\n",Nobs4);
 printf("Total number of Mn2+ in the octahedral lattice: %d\n",Nobs2);
 printf("This corresponds to a fraction of: %lf\n",Nobs2*1.0/((counttype[2])*1.0));

 printf("Leaving replace Mn3 by Mn2 Mn4\n"); 
}


/* Function to exchange some of the Li and Mn */
void swapLitetMnoct(double invfrac)
{
 int i,j,iLi,l1,l2,n,Nobs,Nwanted,dir,NMn,NLidia_beg,NLidia;
 int NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg;
 double randnb;
 FILE *out;

 srand(seed+310);

 /* Number of Mn we want in the tet sites */
 Nwanted=invfrac*(endtet-starttet+1);

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the swap Li-Mn with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the swap Li-Mn with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the swap Li-Mn with NLioctMn3 = %d\n",NLioctMn3_beg); 

 n=0;
 l1=1;
 out=fopen("invfrac_time.dat","w");
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* "Randomly choose iLi */
     iLi=rand()%(endtet-starttet+1)+starttet;
     /* In tet site, if we are not on a Li site, we move to a Li site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[1][iLi]!=1)&&(l2<=lim))
          {
           iLi=rand()%(endtet-starttet+1)+starttet;
           l2++;
          }
     /* Try to jump to one of the neighbouring octahedral sites */
     dir=rand()%12+1;
     /* Swap Litet with Mn2oct */
     if(((bboxtype[2][neighbours_aroundtet[iLi][dir]])==0)&&(neighbours_aroundtet[iLi][dir]<=endoct)&&(xMn<=1.0))
       {
	bboxtype[1][iLi]=bboxtype[2][neighbours_aroundtet[iLi][dir]];
	bboxtype[2][neighbours_aroundtet[iLi][dir]]=1;
	n++;
       }
     /* Swap Litet with Mnoct */
     if(((bboxtype[2][neighbours_aroundtet[iLi][dir]])<=0)&&(neighbours_aroundtet[iLi][dir]<=endoct)&&(xMn>1.0))
       {
	bboxtype[1][iLi]=bboxtype[2][neighbours_aroundtet[iLi][dir]];
	bboxtype[2][neighbours_aroundtet[iLi][dir]]=1;
	n++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
       {
	bboxtype[2][neighbours_aroundtet[iLi][dir]]=bboxtype[1][iLi];
	bboxtype[1][iLi]=1;
	n--;
       }
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
     /*printf("NLidia = %d\n",NLidia); */
     fprintf(out,"%d	%d\n",l1,n);
    }
 fclose(out);
 
 /* Count the actual number of Mn in the Li/tetrahedral lattice */
 Nobs=0;
 for(i=starttet;i<=endtet-starttet+1;i++)
    {
     if(bboxtype[1][i]==0) Nobs++;
    }
 printf("Total number of Mn in the Li/tet lattice: %d\n",Nobs);
 printf("This corresponds to a fraction of: %lf\n",Nobs*1.0/((endtet-starttet+1)*1.0));

}


/* Function to exchange some of the Li and Ti */
void swapLioctTioct(int Nwanted)
{
 int i,j,iLi,iTi,l1,l2,n,Nobs,NMn,NLidia_beg,NLidia;
 int NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg,temp;
 double randnb;

 srand(seed*2);

 /* Number of Mn we want in the tet sites */
 /*Nwanted=invfrac*counttype[1];*/

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the swap Li-Ti in oct sites with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the swap Li-Ti in oct sites with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the swap Li-Ti in oct sites with NLioctMn3 = %d\n",NLioctMn3_beg); 

 /* We don't aim at a specific quantity for now, we just do a certain number of attemps */
 n=1;
 l1=1;
 Nobs=0;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* "Randomly choose iLi */
     iLi=rand()%(endoct-startoct+1)+startoct;
     /* If we are not on a Li site, we move to a Li site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iLi]!=1)&&(l2<=lim))
          {
     	   iLi=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* "Randomly choose iTi */
     iTi=rand()%(endoct-startoct+1)+startoct;
     /* If we are not on a Ti site, we move to a Ti site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iTi]!=2)&&(l2<=lim))
          {
      	   iTi=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Swap Lioct with Tioct */
     if(bboxtype[2][iLi]!=bboxtype[2][iTi])
       {
	temp=bboxtype[2][iLi];
	bboxtype[2][iLi]=bboxtype[2][iTi];
	bboxtype[2][iTi]=temp;
	Nobs++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
       {
	temp=bboxtype[2][iLi];
	bboxtype[2][iLi]=bboxtype[2][iTi];
	bboxtype[2][iTi]=temp;
	Nobs--;
       }
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
     n++;
    }
 
 printf("Number of attemps: %d\n",n);
 printf("Percentage of successful attemps: %lf\n",Nobs*100.0/(n*1.0));

}


/* Function to exchange some of the Li and Mn */
void shuffle_oct(int Nwanted)
{
 int i,j,i1,i2,l1,l2,n,Nobs,NMn,NLidia_beg,NLidia;
 int NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg,temp;
 double randnb;

 srand(seed+49);

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the shuffle Li-Mn in oct sites with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having Mnoct neighbours */
 NLineg_beg=count_Lineg();
 printf("We are starting the shuffle Li-Mn in oct sites with NLineg = %d\n",NLineg_beg); 
 
 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("We are starting the shuffle Li-Mn in oct sites with NLioctMn3 = %d\n",NLioctMn3_beg); 

 /* We don't aim at a specific quantity for now, we just do a certain number of attemps */
 n=1;
 l1=1;
 Nobs=0;
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* "Randomly choose i1 */
     i1=rand()%(endoct-startoct+1)+startoct;
     /* "Randomly choose i2 */
     i2=rand()%(endoct-startoct+1)+startoct;
     /* If we have same type sites, we change. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][i1]==bboxtype[2][i2])&&(l2<=lim))
          {
     	   i2=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Swap the two sites */
     if(bboxtype[2][i1]!=bboxtype[2][i2])
       {
	temp=bboxtype[2][i1];
	bboxtype[2][i1]=bboxtype[2][i2];
	bboxtype[2][i2]=temp;
	Nobs++;
       }
     /* Checking the constraints */
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of NLineg */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
       {
	temp=bboxtype[2][i1];
	bboxtype[2][i1]=bboxtype[2][i2];
	bboxtype[2][i2]=temp;
	Nobs--;
       }
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
     /*printf("NLidia = %d\n",NLidia); */
     n++;
    }
 
 printf("Number of attemps: %d\n",n);
 printf("Percentage of successful attemps: %lf\n",Nobs*1.0/(n*1.0));

}


/* Function to exchange some of the Mntet and Tioct */
void swapLitetTioct(double invfracTi)
{
 int i,j,l1,l2,iLi,iTi,n,Nobs,Nwanted,NLidia,NLidia_beg,dir;
 int NLineg,NLineg_beg,NLioctMn3,NLioctMn3_beg;
 double randnb;
 FILE *out;

 srand(seed+57);

 /* Number of Ti we want in the tet system */
 Nwanted=invfracTi*(endtet-starttet+1);

 /* Count the number of diamagnetic Li at the beginning */
 NLidia_beg=count_Lidia();
 printf("We are starting the swap Litet-Tioct with NLidia = %d\n",NLidia_beg); 

 /* Count the number of Lioct having some Mn in oct */
 NLineg_beg=count_Lineg();
 printf("We are starting the swap Litet-Tioct with NLidia = %d\n",NLineg_beg); 

 /* Count the number of Lioct having Mn3 in oct */
 NLioctMn3_beg=count_LioctMn3();
 printf("Starting with %d Lioct having a Mn3+ in its neighbours\n",NLioctMn3_beg);

 n=0;
 l1=1;
 out=fopen("Tifrac_time.dat","w");
 while((n<Nwanted)&&(l1<=lim))
    {
     l1++;
     /* "Randomly choose iLi */
     iLi=rand()%(endtet-starttet+1)+starttet;
     /* In tet site, if we are not on a Li site, we move to a Li site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[1][iLi]!=1)&&(l2<=lim))
          {
     	   iLi=rand()%(endtet-starttet+1)+starttet;
           l2++;
          }
     /* Randomly choose iTi */
     iTi=rand()%(endoct-startoct+1)+startoct;
     /* In tet site, if we are not on a Li site, we move to a Li site. */
     /* The lim variable serves as a test to not get stuck in the while loop. */
     l2=1;
     while((bboxtype[2][iTi]!=2)&&(l2<=lim))
          {
     	   iTi=rand()%(endoct-startoct+1)+startoct;
           l2++;
          }
     /* Try to jump to one of the neighbouring octahedral sites */
     /*dir=rand()%12+1;*/
     /* Swap Litet with Tioct */
     if(((bboxtype[2][iTi])==2)&&(bboxtype[1][iLi])==1)
       {
	bboxtype[2][iTi]=1;
	bboxtype[1][iLi]=2;
	n++;
       }
     /* Recount the number of diamagnetic Li */
     NLidia=count_Lidia();
     /* Recount the number of Lioct having some Mn in oct */
     NLineg=count_Lineg();
     /* Recount the number of Lioct having Mn3+ in its neighbours */
     NLioctMn3=count_LioctMn3();
     /* If we are going against the constraints, replace back or accept with a given acceptance probability. */
     /* This is done only under constrains */
     randnb=rand()%100000;
     randnb/=100000.0;			/* Get a number between 0 and 1 */
     if(((NLidia>NLidia_beg)||(NLineg>NLineg_beg)||(NLioctMn3>NLioctMn3_beg))&&(randnb>Pacc))
       {
	bboxtype[2][iTi]=2;
	bboxtype[1][iLi]=1;
	n--;
       }
     if(randnb<Pacc)
       {
	NLidia_beg=NLidia;
	NLineg_beg=NLineg;
	NLioctMn3_beg=NLioctMn3;
       }
     /*printf("NLidia = %d\n",NLidia); */
     fprintf(out,"%d	%d\n",l1,n);
    }
 fclose(out);

 /* Count the actual number of Mn in the Ti/tetrahedral lattice */
 Nobs=0;
 for(i=starttet;i<=endtet-starttet+1;i++)
    {
     if(bboxtype[1][i]==2) Nobs++;
    }
 printf("Total number of Ti in the Li/tet lattice: %d\n",Nobs);
 printf("This corresponds to a fraction of: %lf\n",Nobs*1.0/((endtet-starttet+1)*1.0));

}


/* Function to calculate the chemical shifts of Li ions */
void chemical_shift()
{
 FILE *out,*out_tet,*out_oct,*out_p4332,*out_fd3m;
 int i,j,k,NLitet,NLioct,ibin;
 int ncoordMn2,ncoordMn3,ncoordMn4,ncoord2Mn2,ncoord2Mn3,ncoord2Mn4;
 int pncoordMn3,pncoordMn4,fncoordMn3,fncoordMn4,pncoord2Mn3,pncoord2Mn4,fncoord2Mn3,fncoord2Mn4;
 int histoshift[nbins+1],histocoord[12+1][12+1],histocoord2[6+1][6+1];
 int histoshift_tet[nbins+1],histoshift_oct[nbins+1];
 int histoshift_p4332[nbins+1],histoshift_fd3m[nbins+1];
 double Lishift,avshift_tet,avshift_oct,dshift,dx,dy,dz,dist,maxint,integral,avoxstate,pavoxstate,favoxstate;
 double pLitet_Mn3oct,pLitet_Mn4oct,fLitet_Mn3oct,fLitet_Mn4oct,pLioct_Mn3oct,pLioct_Mn4oct,fLioct_Mn3oct,fLioct_Mn4oct;
 double *modelNMR,*modelNMR_bis,*modelNMR_tet,*modelNMR_oct,*modelNMR_p4332,*modelNMR_fd3m;
 double ***modelNMR_tet_decomp,***modelNMR_oct_decomp;

 /* IMPORTANT */
 /* For xMn<=1.0, P4332 and FD3M are P4332 and FD3M
    For xMn>1.0, P4332 is the Ti-rich region and FD3M is the Mn-rich region */

 /* Shifts (ppm) obtained from HYB20 calculations */
 if(strcmp(shifttype,"HYB20")==0) 
   {Lioct_Mn2tet=18.1; Lioct_Mn2oct=-72.5; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4; 
    Litet_Mn2oct=31.6; Litet_Mn3oct=31.3; Litet_Mn4oct=63.72;} 

 /* Shifts (ppm) obtained from HYB35 calculations */
 if(strcmp(shifttype,"HYB35")==0) 
   {Lioct_Mn2tet=15.3; Lioct_Mn2oct=-54.1; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4; 
    Litet_Mn2oct=23.7; Litet_Mn3oct=26.2; Litet_Mn4oct=50.25;} 

 /* Shifts (ppm) obtained averages and modifications */
 if(strcmp(shifttype,"AVDFT")==0) 
   {Lioct_Mn2tet=16.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;
    Litet_Mn2oct=27.65; Litet_Mn3oct=28.8; Litet_Mn4oct=56.985;} 
 if(strcmp(shifttype,"AVmodif")==0) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;
    Litet_Mn2oct=32.65; Litet_Mn3oct=28.8; Litet_Mn4oct=56.985;} 
 /* Shifts (ppm) obtained for dynamic averaging of Mn3+ and Mn4+ in oct */
 /* AVdynam are taken following AVmodif */
 /* First case, there is actually only one phase */
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac==0)||(ffrac==0))) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;
    Litet_Mn2oct=32.65; Litet_Mn3oct=28.8; Litet_Mn4oct=56.985;
    avoxstate=(4.0*xMn-1.0)/xMn; 
    /* if oxredfrac!=0, need to update avoxstate (reg1 or reg2 will be empty) */
    if(oxredfrac>0) avoxstate=(4.0*(reg1_NMn4+reg2_NMn4)+3.0*(reg1_NMn3+reg2_NMn3))/(reg1_NMn4+reg2_NMn4+reg1_NMn3+reg2_NMn3);
    Litet_Mn3oct=(4.0-avoxstate)*Litet_Mn3oct+(avoxstate-3.0)*Litet_Mn4oct; Litet_Mn4oct=Litet_Mn3oct;
    Lioct_Mn3oct=(4.0-avoxstate)*Lioct_Mn3oct+(avoxstate-3.0)*Lioct_Mn4oct; Lioct_Mn4oct=Lioct_Mn3oct;}
 /* Second case, there are two phases */
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac!=0)&&(ffrac!=0))) 
   {Lioct_Mn2tet=21.7; Lioct_Mn2oct=-63.3; Lioct_Mn3oct=124.1; Lioct_Mn4oct=181.4;
    Litet_Mn2oct=32.65; Litet_Mn3oct=28.8; Litet_Mn4oct=56.985;
    pavoxstate=(4.0*pxMn-1.0)/pxMn; 
    favoxstate=(4.0*fxMn-1.0)/fxMn; 
    /* if oxredfrac!=0, need to update avoxstate */
    if(oxredfrac>0) pavoxstate=(4.0*reg1_NMn4+3.0*reg1_NMn3)/(reg1_NMn4+reg1_NMn3);
    if(oxredfrac>0) favoxstate=(4.0*reg2_NMn4+3.0*reg2_NMn3)/(reg2_NMn4+reg2_NMn3);
    pLitet_Mn3oct=(4.0-pavoxstate)*Litet_Mn3oct+(pavoxstate-3.0)*Litet_Mn4oct; pLitet_Mn4oct=pLitet_Mn3oct;
    pLioct_Mn3oct=(4.0-pavoxstate)*Lioct_Mn3oct+(pavoxstate-3.0)*Lioct_Mn4oct; pLioct_Mn4oct=pLioct_Mn3oct;
    fLitet_Mn3oct=(4.0-favoxstate)*Litet_Mn3oct+(favoxstate-3.0)*Litet_Mn4oct; fLitet_Mn4oct=fLitet_Mn3oct;
    fLioct_Mn3oct=(4.0-favoxstate)*Lioct_Mn3oct+(favoxstate-3.0)*Lioct_Mn4oct; fLioct_Mn4oct=fLioct_Mn3oct;}
 
 /* Set histoshift to zero */
 for(i=0;i<=nbins;i++) {histoshift[i]=0; histoshift_tet[i]=0; histoshift_oct[i]=0; histoshift_p4332[i]=0; histoshift_fd3m[i]=0;}
 modelNMR=dvector(nbins2+1);
 modelNMR_bis=dvector(nbins2+1);
 modelNMR_tet=dvector(nbins2+1);
 modelNMR_oct=dvector(nbins2+1);
 modelNMR_p4332=dvector(nbins2+1);
 modelNMR_fd3m=dvector(nbins2+1);
 modelNMR_tet_decomp=tddmatrix(12+1,12+1,nbins2+1);
 modelNMR_oct_decomp=tddmatrix(6+1,6+1,nbins2+1);
 for(i=0;i<=nbins2;i++) 
    {
     modelNMR[i]=0.0; modelNMR_bis[i]=0.0; modelNMR_tet[i]=0.0; modelNMR_oct[i]=0.0; modelNMR_p4332[i]=0.0; modelNMR_fd3m[i]=0.0;
     for(j=0;j<=12;j++)
	{
	 for(k=0;k<=12;k++) 
	    {
	     modelNMR_tet_decomp[j][k][i]=0.0; 
	     if((j<=6)&&(k<=6)) modelNMR_oct_decomp[j][k][i]=0.0;
	    }
	}
    }

 dshift=(smax-smin)/(nbins*1.0);

 /* Only one phase or not using dynamical averaging for Mn3+/Mn4+ */
 if((strcmp(shifttype,"AVdynam")!=0)||((pfrac==0)||(ffrac==0))) 
   {
    /* First check the Li in tet sites */
    avshift_tet=0.0;
    NLitet=0;
    for(i=0;i<=12;i++) {for(j=0;j<=12;j++) histocoord[i][j]=0;}
    for(i=1;i<=counttype[1];i++)
    {
     /* If there is a lithium, check the number of Mn in oct */
     if(bboxtype[1][i]==1)
       {
	NLitet++;
	ncoordMn2=0; ncoordMn3=0; ncoordMn4=0;
        for(j=1;j<=counttype[2];j++)
	   {
	    dx=posx[1][i]-posx[2][j]; dy=posy[1][i]-posy[2][j]; dz=posz[1][i]-posz[2][j];
            /* Periodic boundary conditions */
            if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
            if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
            if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	    dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	    if((dist<=Rcut)&&(bboxtype[2][j]==0)) {ncoordMn2++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-1)) {ncoordMn3++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-2)) {ncoordMn4++;}
	   }
	if(xMn<=1.0) histocoord[ncoordMn2][ncoordMn3]++;
	if(xMn>1.0) histocoord[ncoordMn4][ncoordMn3]++;
	avshift_tet+=ncoordMn2*Litet_Mn2oct+ncoordMn3*Litet_Mn3oct+ncoordMn4*Litet_Mn4oct;
	Lishift=ncoordMn2*Litet_Mn2oct+ncoordMn3*Litet_Mn3oct+ncoordMn4*Litet_Mn4oct;
        ibin=floor((Lishift-smin)/dshift);
        if((ibin>=0)&&(ibin<=nbins))
          {
           histoshift[ibin]++;
           histoshift_tet[ibin]++;
	   if(i<=(counttype[1]*pfrac)) histoshift_p4332[ibin]++;
	   if(i>(counttype[1]*pfrac)) histoshift_fd3m[ibin]++;
          }
       }
    }

    /* Then check the Li in oct sites */
    avshift_oct=0.0;
    NLioct=0;
    for(i=0;i<=6;i++) 
    {
     for(j=1;j<=6;j++) histocoord2[i][j]=0; 
    }
    for(i=1;i<=counttype[2];i++)
    {
     /* If there is a lithium, check the number of Mn in tet and in oct */
     if(bboxtype[2][i]==1)
       {
	NLioct++;
	/* Mn in tetrahedral site, only Mn2+ is allowed to go in tet site */
	ncoordMn2=0;
        for(j=1;j<=counttype[1];j++)
	   {
	    dx=posx[2][i]-posx[1][j]; dy=posy[2][i]-posy[1][j]; dz=posz[2][i]-posz[1][j]; 
	    /* Periodic boundary conditions */
            if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
            if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
            if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	    dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	    if((dist<=Rcut)&&(bboxtype[1][j]==0)) {ncoordMn2++;}
	   }
	avshift_oct+=ncoordMn2*Lioct_Mn2tet;	
	Lishift=ncoordMn2*Lioct_Mn2tet;
	/* Mn in octahedral site */
	ncoord2Mn2=0; ncoord2Mn3=0; ncoord2Mn4=0;
        for(j=1;j<=counttype[2];j++)
	   {
	    if(j!=i)
	      {
	       dx=posx[2][i]-posx[2][j]; dy=posy[2][i]-posy[2][j]; dz=posz[2][i]-posz[2][j];
               /* Periodic boundary conditions */
               if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
               if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
               if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	       dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	       if((dist<=Rcut)&&(bboxtype[2][j]==0)) {ncoord2Mn2++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-1)) {ncoord2Mn3++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-2)) {ncoord2Mn4++;}
	      }
	   }
	if(xMn<=1.0) histocoord2[ncoordMn2][ncoord2Mn2]++;
	if(xMn>1.0) histocoord2[ncoord2Mn3][ncoord2Mn4]++;
	avshift_oct+=ncoord2Mn2*Lioct_Mn2oct+ncoord2Mn3*Lioct_Mn3oct+ncoord2Mn4*Lioct_Mn4oct;	
	Lishift+=ncoord2Mn2*Lioct_Mn2oct+ncoord2Mn3*Lioct_Mn3oct+ncoord2Mn4*Lioct_Mn4oct;
	/* Finally do the histogram of shifts */
        ibin=floor((Lishift-smin)/dshift);
        if((ibin>=0)&&(ibin<=nbins))
          {
           histoshift[ibin]++;
           histoshift_oct[ibin]++;
	   if(i<=(counttype[2]*pfrac)) histoshift_p4332[ibin]++;
	   if(i>(counttype[2]*pfrac)) histoshift_fd3m[ibin]++;
          }
       }
    }
   }

 /* Two phases and using dynamical averaging for Mn3+/Mn4+ */
 if((strcmp(shifttype,"AVdynam")==0)&&((pfrac!=0)&&(ffrac!=0))) 
   {
    /* First check the Li in tet sites */
    avshift_tet=0.0;
    NLitet=0;
    for(i=0;i<=12;i++) {for(j=0;j<=12;j++) histocoord[i][j]=0;}
    for(i=1;i<=counttype[1];i++)
    {
     /* If there is a lithium, check the number of Mn in oct */
     if(bboxtype[1][i]==1)
       {
	NLitet++;
	ncoordMn2=0; ncoordMn3=0; ncoordMn4=0;
	pncoordMn3=0; pncoordMn4=0;
	fncoordMn3=0; fncoordMn4=0;
        for(j=1;j<=counttype[2];j++)
	   {
	    dx=posx[1][i]-posx[2][j]; dy=posy[1][i]-posy[2][j]; dz=posz[1][i]-posz[2][j];
            /* Periodic boundary conditions */
            if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
            if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
            if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	    dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	    if((dist<=Rcut)&&(bboxtype[2][j]==0)) {ncoordMn2++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-1)&&(j<=counttype[2]*pfrac)) {ncoordMn3++; pncoordMn3++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-2)&&(j<=counttype[2]*pfrac)) {ncoordMn4++; pncoordMn4++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-1)&&(j>counttype[2]*pfrac)) {ncoordMn3++; fncoordMn3++;}
	    if((dist<=Rcut)&&(bboxtype[2][j]==-2)&&(j>counttype[2]*pfrac)) {ncoordMn4++; fncoordMn4++;}
	   }
	if(xMn<=1.0) histocoord[ncoordMn2][ncoordMn3]++;
	if(xMn>1.0) histocoord[ncoordMn4][ncoordMn3]++;
	avshift_tet+=ncoordMn2*Litet_Mn2oct+pncoordMn3*pLitet_Mn3oct+pncoordMn4*pLitet_Mn4oct+fncoordMn3*fLitet_Mn3oct+fncoordMn4*fLitet_Mn4oct;
	Lishift=ncoordMn2*Litet_Mn2oct+pncoordMn3*pLitet_Mn3oct+pncoordMn4*pLitet_Mn4oct+fncoordMn3*fLitet_Mn3oct+fncoordMn4*fLitet_Mn4oct;
        ibin=floor((Lishift-smin)/dshift);
        if((ibin>=0)&&(ibin<=nbins))
          {
           histoshift[ibin]++;
           histoshift_tet[ibin]++;
	   if(i<=(counttype[1]*pfrac)) histoshift_p4332[ibin]++;
	   if(i>(counttype[1]*pfrac)) histoshift_fd3m[ibin]++;
          }
       }
    }

    /* Then check the Li in oct sites */
    avshift_oct=0.0;
    NLioct=0;
    for(i=0;i<=6;i++) 
    {
     for(j=1;j<=6;j++) histocoord2[i][j]=0; 
    }
    for(i=1;i<=counttype[2];i++)
    {
     /* If there is a lithium, check the number of Mn in tet and in oct */
     if(bboxtype[2][i]==1)
       {
	NLioct++;
	/* Mn in tetrahedral site, only Mn2+ is allowed to go in tet site */
	ncoordMn2=0;
        for(j=1;j<=counttype[1];j++)
	   {
	    dx=posx[2][i]-posx[1][j]; dy=posy[2][i]-posy[1][j]; dz=posz[2][i]-posz[1][j]; 
	    /* Periodic boundary conditions */
            if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
            if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
            if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	    dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	    if((dist<=Rcut)&&(bboxtype[1][j]==0)) {ncoordMn2++;}
	   }
	avshift_oct+=ncoordMn2*Lioct_Mn2tet;	
	Lishift=ncoordMn2*Lioct_Mn2tet;
	/* Mn in octahedral site */
	ncoord2Mn2=0; ncoord2Mn3=0; ncoord2Mn4=0;
	pncoord2Mn3=0; pncoord2Mn4=0; fncoord2Mn3=0; fncoord2Mn4=0;
        for(j=1;j<=counttype[2];j++)
	   {
	    if(j!=i)
	      {
	       dx=posx[2][i]-posx[2][j]; dy=posy[2][i]-posy[2][j]; dz=posz[2][i]-posz[2][j];
               /* Periodic boundary conditions */
               if(dx>(Lx/2.0)) dx-=Lx; if(dx<(-Lx/2.0)) dx+=Lx;
               if(dy>(Ly/2.0)) dy-=Ly; if(dy<(-Ly/2.0)) dy+=Ly;
               if(dz>(Lz/2.0)) dz-=Lz; if(dz<(-Lz/2.0)) dz+=Lz;
	       dist=dx*dx+dy*dy+dz*dz; dist=sqrt(dist);
	       if((dist<=Rcut)&&(bboxtype[2][j]==0)) {ncoord2Mn2++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-1)&&(j<=counttype[2]*pfrac)) {ncoord2Mn3++; pncoord2Mn3++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-2)&&(j<=counttype[2]*pfrac)) {ncoord2Mn4++; pncoord2Mn4++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-1)&&(j>counttype[2]*pfrac)) {ncoord2Mn3++; fncoord2Mn3++;}
	       if((dist<=Rcut)&&(bboxtype[2][j]==-2)&&(j>counttype[2]*pfrac)) {ncoord2Mn4++; fncoord2Mn4++;}
	      }
	   }
	if(xMn<=1.0) histocoord2[ncoordMn2][ncoord2Mn2]++;
	if(xMn>1.0) histocoord2[ncoord2Mn3][ncoord2Mn4]++;
	avshift_oct+=ncoord2Mn2*Lioct_Mn2oct+pncoord2Mn3*pLioct_Mn3oct+pncoord2Mn4*pLioct_Mn4oct+fncoord2Mn3*fLioct_Mn3oct+fncoord2Mn4*fLioct_Mn4oct;	
	Lishift+=ncoord2Mn2*Lioct_Mn2oct+pncoord2Mn3*pLioct_Mn3oct+pncoord2Mn4*pLioct_Mn4oct+fncoord2Mn3*fLioct_Mn3oct+fncoord2Mn4*fLioct_Mn4oct;	
	/* Finally do the histogram of shifts */
        ibin=floor((Lishift-smin)/dshift);
        if((ibin>=0)&&(ibin<=nbins))
          {
           histoshift[ibin]++;
           histoshift_oct[ibin]++;
	   if(i<=(counttype[2]*pfrac)) histoshift_p4332[ibin]++;
	   if(i>(counttype[2]*pfrac)) histoshift_fd3m[ibin]++;
          }
       }
    }
   }

 /* Print the coordination histograms and calculate the "decomposed" model NMR spectra */ 
 printf("\nFor Li in tetrahedral sites:\n");
 if(NLitet!=0)
   {
    printf("Average shift = %lf ppm\n",avshift_tet/(NLitet*1.0));
    out=fopen("Histocoord_Litet.dat","w");
    if(xMn<=1.0) fprintf(out,"# Nb of Mn2+ in the coordination shell - Nb of Mn3+ in the coordination shell - percentage of Li with this number - shift (ppm)\n");
    if(xMn>1.0) fprintf(out,"# Nb of Mn4+ in the oct coordination shell - Nb of Mn3+ in the coordination shell - percentage of Li with this number - shift (ppm)\n");
    for(i=0;i<=12;i++) 
       {
	for(j=0;j<=12;j++)
	   {
	    if(xMn<=1.0) Lishift=i*Litet_Mn2oct+j*Litet_Mn3oct;
	    if(xMn>1.0) Lishift=i*Litet_Mn4oct+j*Litet_Mn3oct;
	    if(histocoord[i][j]!=0) 
	      {
	       fprintf(out,"%d	%d	%lf	%lf\n",i,j,histocoord[i][j]*100.0/((NLitet+NLioct)*1.0),Lishift);
     	       for(k=0;k<=nbins2;k++)
        	  {
	 	   modelNMR_tet_decomp[i][j][k]+=histocoord[i][j]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(k*(smax-smin)/(nbins2*1.0)+smin-(Lishift))*(k*(smax-smin)/(nbins2*1.0)+smin-(Lishift))/(2*gwidth*gwidth));
        	  }
	      }
	   }
       }
    fclose(out);
   }
 if(NLitet==0)
   {
    printf("No Li in tetrahedral sites\n");
   }

 printf("\nFor Li in octahedral sites:\n");
 if(NLioct!=0)
   {
    printf("Average shift = %lf ppm\n",avshift_oct/(NLioct*1.0));
    out=fopen("Histocoord_Lioct.dat","w");
    if(xMn<=1.0) fprintf(out,"# Number of Mn2 in the tet shell - Number of Mn2 in the oct shell - percentage of Li having this configuration - shift (ppm)\n");
    if(xMn>1.0) fprintf(out,"# Number of Mn3 in the oct shell - Number of Mn4 in the oct shell - percentage of Li having this configuration - shift (ppm)\n");
    for(i=0;i<=6;i++) 
       {
        for(j=0;j<=6;j++) 
	   {
	    Lishift=i*Lioct_Mn2tet+j*Lioct_Mn2oct;
	    if(histocoord2[i][j]!=0) 
	      {
	       fprintf(out,"%d	%d	%lf	%lf\n",i,j,histocoord2[i][j]*100.0/((NLitet+NLioct*1.0)),Lishift);
     	       for(k=0;k<=nbins2;k++)
        	  {
	 	   modelNMR_oct_decomp[i][j][k]+=histocoord2[i][j]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(k*(smax-smin)/(nbins2*1.0)+smin-(Lishift))*(k*(smax-smin)/(nbins2*1.0)+smin-(Lishift))/(2*gwidth*gwidth));
        	  }
	      }
	   }
       }
    fclose(out);
   }
 if(NLioct==0)
   {
    printf("No Li in octahedral sites\n");
   }

 /* Write down the histogram of shifts + calculate the model NMR spectrum */
 out=fopen("Histo_shift.dat","w");
 out_tet=fopen("Histo_shift_tet.dat","w");
 out_oct=fopen("Histo_shift_oct.dat","w");
 if(xMn<=1.0) 
   {out_p4332=fopen("Histo_shift_P4332.dat","w");
    out_fd3m=fopen("Histo_shift_Fd3m.dat","w");}
 if(xMn>1.0) 
   {out_p4332=fopen("Histo_shift_Ti-rich.dat","w");
    out_fd3m=fopen("Histo_shift_Mn-rich.dat","w");}
 fprintf(out,"# Shift (ppm) - percentage of Li with this shift\n");
 fprintf(out_tet,"# Shift (ppm) - percentage of Li with this shift\n");
 fprintf(out_oct,"# Shift (ppm) - percentage of Li with this shift\n");
 fprintf(out_p4332,"# Shift (ppm) - percentage of Li with this shift\n");
 fprintf(out_fd3m,"# Shift (ppm) - percentage of Li with this shift\n");
 for(i=0;i<=nbins;i++)
    {
     if(histoshift[i]!=0) fprintf(out,"%lf      %lf\n",i*dshift+smin,histoshift[i]*100.0/((NLitet+NLioct)*1.0));
     if(histoshift_tet[i]!=0) fprintf(out_tet,"%lf      %lf\n",i*dshift+smin,histoshift_tet[i]*100.0/((NLitet+NLioct)*1.0));
     if(histoshift_oct[i]!=0) fprintf(out_oct,"%lf      %lf\n",i*dshift+smin,histoshift_oct[i]*100.0/((NLitet+NLioct)*1.0));
     if(histoshift_p4332[i]!=0) fprintf(out_p4332,"%lf      %lf\n",i*dshift+smin,histoshift_p4332[i]*100.0/((NLitet+NLioct)*1.0));
     if(histoshift_fd3m[i]!=0) fprintf(out_fd3m,"%lf      %lf\n",i*dshift+smin,histoshift_fd3m[i]*100.0/((NLitet+NLioct)*1.0));
     for(j=0;j<=nbins2;j++)
        {
	 modelNMR[j]+=histoshift_p4332[i]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth*gwidth));
	 modelNMR[j]+=histoshift_fd3m[i]*100.0/((NLitet+NLioct)*1.0*gwidth2)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth2*gwidth2));
	 modelNMR_bis[j]+=histoshift[i]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth*gwidth));
	 modelNMR_tet[j]+=histoshift_tet[i]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth*gwidth));
	 modelNMR_oct[j]+=histoshift_oct[i]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth*gwidth));
	 modelNMR_p4332[j]+=histoshift_p4332[i]*100.0/((NLitet+NLioct)*1.0*gwidth)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth*gwidth));
	 modelNMR_fd3m[j]+=histoshift_fd3m[i]*100.0/((NLitet+NLioct)*1.0*gwidth2)*exp(-1.0*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))*(j*(smax-smin)/(nbins2*1.0)+smin-(i*dshift+smin))/(2*gwidth2*gwidth2));
	}
    }
 fclose(out);
 fclose(out_tet);
 fclose(out_oct);
 fclose(out_p4332);
 fclose(out_fd3m);

 printf("\nNLitot / NLitet / NLioct %d %d %d\n",counttype[1],NLitet,NLioct);

 /* Find the maximum intensity for the normalisation of the spectrum */
 maxint=0.0;
 integral=0.0;
 for(j=0;j<=nbins2;j++)
    {
     if(j<nbins2) integral+=(modelNMR[j]+modelNMR[j+1])*0.5*dshift2;
     if(maxint<=modelNMR[j]) maxint=modelNMR[j]; 
    }
 printf("The integral is: %lf\n",integral);

 /* Write down the model NMR spectrum */
 out=fopen("Model_NMR_spectrum.dat","w");
 out_tet=fopen("Model_NMR_spectrum_tet.dat","w");
 out_oct=fopen("Model_NMR_spectrum_oct.dat","w");
 if(xMn<=1.0)
   {out_p4332=fopen("Model_NMR_spectrum_P4332.dat","w");
    out_fd3m=fopen("Model_NMR_spectrum_Fd3m.dat","w");}
 if(xMn>1.0)
   {out_p4332=fopen("Model_NMR_spectrum_Ti-rich.dat","w");
    out_fd3m=fopen("Model_NMR_spectrum_Mn-rich.dat","w");}
 fprintf(out,"# Shift (ppm) - Intensity\n");
 fprintf(out_tet,"# Shift (ppm) - Intensity\n");
 fprintf(out_oct,"# Shift (ppm) - Intensity\n");
 fprintf(out_p4332,"# Shift (ppm) - Intensity\n");
 fprintf(out_fd3m,"# Shift (ppm) - Intensity\n");
 for(j=0;j<=nbins2;j++)
    { 
     fprintf(out,"%lf	%lf	%lf\n",j*(smax-smin)/(nbins2*1.0)+smin,modelNMR[j]/maxint,modelNMR_bis[j]/maxint);
     fprintf(out_tet,"%lf	%lf\n",j*(smax-smin)/(nbins2*1.0)+smin,modelNMR_tet[j]/maxint);
     fprintf(out_oct,"%lf	%lf\n",j*(smax-smin)/(nbins2*1.0)+smin,modelNMR_oct[j]/maxint);
     fprintf(out_p4332,"%lf	%lf\n",j*(smax-smin)/(nbins2*1.0)+smin,modelNMR_p4332[j]/maxint);
     fprintf(out_fd3m,"%lf	%lf\n",j*(smax-smin)/(nbins2*1.0)+smin,modelNMR_fd3m[j]/maxint);
    }
 fclose(out);
 fclose(out_tet);
 fclose(out_oct);
 fclose(out_p4332);
 fclose(out_fd3m);

 /* Write down the decomposed spectra */
 out_tet=fopen("Decomp_model_NMR_spectrum_tet.dat","w");
 out_oct=fopen("Decomp_model_NMR_spectrum_oct.dat","w");
 fprintf(out_tet,"# Shift (ppm) - Intensity\n");
 fprintf(out_oct,"# Shift (ppm) - Intensity\n");
 for(k=0;k<=nbins2;k++)
    {
     fprintf(out_tet,"%4.4lf	",k*(smax-smin)/(nbins2*1.0)+smin);
     fprintf(out_oct,"%4.4lf	",k*(smax-smin)/(nbins2*1.0)+smin);
     for(i=0;i<=12;i++)
        { 
	 for(j=0;j<=12;j++)
	    {
	     if(histocoord[i][j]!=0) fprintf(out_tet,"%4.4lf	",modelNMR_tet_decomp[i][j][k]/maxint);
	     if((histocoord2[i][j]!=0)&&(i<=6)&&(j<=6)) fprintf(out_oct,"%4.4lf	",modelNMR_oct_decomp[i][j][k]/maxint);
	    }
        }
     fprintf(out_tet,"\n");
     fprintf(out_oct,"\n");
    }
 fclose(out_tet);
 fclose(out_oct);

}


/******************************* Dynamic memory allocation *****************************/

int **imatrix(int nl, int nc)
{
  int i;
  int **m;
  m=(int **) malloc(nl*sizeof(int*));
  if (m) { m[0]=(int *) malloc(nl*nc*sizeof(int));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

double *dvector(int n)
{
  return (double *) malloc(n*sizeof(double));
}

double **dmatrix(int nl, int nc)
{
  int i;
  double **m;
  m=(double **) malloc(nl*sizeof(double*));
  if (m) { m[0]=(double *) malloc(nl*nc*sizeof(double));
           if (m[0]==NULL) return NULL;
           for (i=1;i<nl;i++) m[i]=m[i-1]+nc;
         }
  return m;
}

int *ivector(int n)
{
  return (int *) malloc(n*sizeof(int));
}

double ***tddmatrix (int X_SIZE, int Y_SIZE, int Z_SIZE)
{
 double ***m;
 int i,j;

 m = (double ***)malloc(sizeof(double **) * X_SIZE);
  
 for (i = 0 ;  i < X_SIZE; i++) 
 {
    m[i] = (double **)malloc(sizeof(double *) * Y_SIZE);
  
    for (j = 0; j < Y_SIZE; j++)
       m[i][j] = (double *)malloc(sizeof(double) * Z_SIZE);
 }

 return m;
}

