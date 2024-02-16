/*Read me*/

/* 
The following code is an example on how to perform 2d numerical simulation of a reaction-diffusion system obeying the
normal form using a 4th-order Runge-Kutta algorithm.
To compile it, you need a .c compiler on your machine, for instance gcc.
In the terminal, type: "gcc Normal_form_2d.c -lm -fopenmp", to compile the code.
Then: "./a.out", to execute it.
Output files are .dat.

-lm: calls the math library of your machine
-fopenmp: used for multithreading (enables pragma)

This codes can uses gasdev or ran1 to generate random numbers for the initial conditions.

*/
/***********************************************************************************************/
//These are the libraries
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
/***********************************************************************************************/
//The following allows to work with big matrices by properly allocating pointers and rows
#define NR_END 1
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
fprintf(stderr,"Numerical Recipes run-time error...\n");
fprintf(stderr,"%s\n",error_text);
fprintf(stderr,"...now exiting to system...\n");
exit(1);
}
double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
double **m;
/* allocate pointers to rows */
m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;
/* allocate rows and set pointers to them */
m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;
for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
/* return pointer to array of pointers to rows */
return m;
}
/***********************************************************************************************/
// These are the files
FILE *fichier2, *fichier3, *fichierCtot;
FILE *fichiersU[1000];
FILE *fichierCI;
char filenameCI[200];
char filename1[200], filename2[200], filenameCtot[200];
//These are the concentration pointers and vectors
double **CU,**CCU,**vt_CU;
double CU_middle,CU_minusi ,CU_plusi ,CU_minusj,CU_plusj;
//Ran1 & gasdev
#define ttot 10000.
long m, ti1, tim, *seed;
float ran1(long *idum);
float gasdev(long *idum);
double r1,g1;
//Some definitions
double resU,K1U,K2U,K3U,K4U;
double Ucount,SumU,old_SumU;
int numsnap,s,snap,tprintCtot;
int N,T,t, tprint, nprint,i,j,iprint,Sumprint,k,pmax,p,tsnap;
double step;                        /*step = D / dx²*/                                 
/*----------------------------------------------< End: declaration >---------------------------------------------------------------------------------------------*/
/*Main*/
int main()
{
  int NumSnap = 1000; //# of snapshots, i.e. files containing the concentration maps
  double L = 100.;  //size of the system
  double dx = 1.; //spatial increment
  double Ttot = 20000.; //total simulation time
  double dt = 0.1;  //time increment
  //Files creation and opening
  for (int s = 0; s < NumSnap; s++)
    {
        char filenameU[200];
        sprintf(filenameU, "Map_U_0%d.dat", s);
        fichiersU[s] = fopen(filenameU, "w");
    }                                 
    sprintf(filenameCtot,"Total_concentration.dat");      
    sprintf(filenameCI,"Seed_for_initial_condition.dat");       
    fichierCtot = fopen(filenameCtot, "w");
    fichierCI = fopen(filenameCI, "w");
    //Number of integration steps and element of the matrix
      N = round(L/dx);
      T = round(Ttot/dt);
      printf("Matrix : %d ; #time step : %d \n",N*N,T);
      //Matrix
      CU = dmatrix(0,N,0,N);
      CCU = dmatrix(0,N,0,N);
      vt_CU = dmatrix(0,N,0,N);
      //Random number generation
      m = -time(&tim);
      seed = &m;
      ti1 = time(&tim) % 549874;
      for (p=0; p < ti1; p++){
          ran1(seed);  //random distribution
      }
      for (j=0; j < ti1; j++){
          gasdev(seed); //gaussian distribution
      }
      //setting to 0
      Ucount = 0;
      tprint = 0;
      tprintCtot = 0;
      snap = 0;
      //step
      step = 1/(dx*dx);
      //parameters
      double lambda,g,gamma,up,un;
      lambda = 0.5;
      g = 0.;
      gamma = 0.;
      up = sqrt(lambda);
      un = -sqrt(lambda);
      /***********************************************************************************************/
      /*Initial conditions*/
      for (i = 1 ; i < N ; i++)
      {
        for (j = 1; j < N; j++)
        {
          g1 = gasdev(seed);
          fprintf(fichierCI,"%lf \n",g1); // --> Here we use gaussian random numbers for the initial condition
          CU[i][j] = g1*0.01;
          Ucount += g1*0.01;
        }
      }
      fclose(fichierCI);
      /***********************************************************************************************/
      /*Integration loops*/
      for(t = 0 ; t <= T ; t++)// Time loop 
      {
        if (t == tprintCtot)
        {
          printf("[%d/%d] \n",tprintCtot,T); //Print the progress of the simulation in time
        }
        /*Boundary conditions*/
          for (i = 0; i <= N; i++)
          {
            CU[i][0] = (4.*CU[i][1] - CU[i][2])/3.;
            CU[i][N] = (4.*CU[i][N-1] - CU[i][N-2])/3.;  
          }
          for (j = 0; j <= N; j++)
          {
            CU[0][j] = (4.*CU[1][j] - CU[2][j])/3.;
            CU[N][j] = (4.*CU[N-1][j] - CU[N-2][j])/3.;
          }
          //The following command 'pragma' allows to take advantage of multithreading to accelerate the simulation computation time
          #pragma omp parallel for num_threads(20) private (Ucount,resU,j,t,iprint,K1U,K2U,K3U,K4U,r1,g1,CU_middle,CU_minusi,CU_plusi,CU_minusj,CU_plusj)
          for (i = 1; i < N; i++) // 2D spatial loops
          {
            for (j = 1; j < N; j++)
            {
              //Assignment vector --> pointer to save computational efforts
              CU_middle = CU[i][j];
              CU_minusi = CU[i-1][j];
              CU_plusi = CU[i+1][j];
              CU_minusj = CU[i][j-1];
              CU_plusj = CU[i][j+1];
              /* Runge-Kutta 4 */
              K1U = dt*(CU_middle*lambda - CU_middle*CU_middle*CU_middle + g);
              K2U = dt*((CU_middle + K1U/2)*lambda - (CU_middle + K1U/2)*(CU_middle + K1U/2)*(CU_middle + K1U/2) + g);
              K3U = dt*((CU_middle + K2U/2)*lambda - (CU_middle + K2U/2)*(CU_middle + K2U/2)*(CU_middle + K2U/2) + g);      
              K4U = dt*((CU_middle + K3U)*lambda - (CU_middle + K3U)*(CU_middle + K3U)*(CU_middle + K3U) + g);
              // This the reaction term
              resU = CU_middle + K1U/6 + K2U/3 + K3U/3 + K4U/6;
              // This is the reaction term + the diffusion contribution
              CCU[i][j] = resU + dt*(((1. + 2.*gamma*CU_middle)*(CU_plusi + CU_minusi - 2.*CU_middle)/(dx*dx)) + gamma*2.*((CU_plusi - CU_minusi)/(2.*dx))*((CU_plusi - CU_minusi)/(2.*dx))) + dt*(((1. + 2.*gamma*CU_middle)*(CU_plusj + CU_minusj - 2.*CU_middle)/(dx*dx)) + gamma*2.*((CU_plusj - CU_minusj)/(2.*dx))*((CU_plusj - CU_minusj)/(2.*dx)));
            }
          }//End of spatial loops 

          /***********************************************************************************************/

          /* Variables updates and results printing*/
          old_SumU = SumU;
          SumU = 0;
          //The following loop calculate the average concentration of u and print a concentration map for a given timestep
          for (i = 1; i < N; i++)
          {
            for (j = 1; j < N; j++)
            {
              SumU += CU[i][j]/((N-2)*(N-2));
              if (t == tprint) 
                  {
                    fprintf(fichiersU[snap], "%lf \t",CU[i][j]);
                  }      
              }
            if (t == tprint)
              { 
                fprintf(fichiersU[snap], "\n");  
              } 
          }
        if (t == tprint)
        {
          tprint += (int)(T/NumSnap);
          fclose(fichiersU[snap]);
          snap += 1;
        }
          //Update of the concentration pointers for the next integration
          vt_CU = CU;
          CU = CCU;
          CCU = vt_CU;
         if (t == tprintCtot)
        {
          fprintf(fichierCtot, "%lf \t %lf \n", (dt*t), SumU);
          tprintCtot += (int)(T/500);
        }
      }  //End of time loop

    return 0;
}

/* End of program */ 

// Ran1 for random numbers

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

// Box-Muller technique for Gaussian random numbers

#include <math.h>
float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;
  if (iset == 0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
    
  }
}


/*------------------------------------------------------------------< End >---------------------------------------------------------------------------------------*/