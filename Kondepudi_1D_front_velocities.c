/*Read me*/

/* 
The following code is an example on how to perform multiple 1d numerical simulation of a reaction-diffusion system obeying the
Kondepudi model to study the propagation velocity of a front separating two domains of opposite chirality.
To compile it, you need a .c compiler on your machine, for instance gcc.
In the terminal, type: "gcc Kondepudi_1D_front_velocities.c -lm -fopenmp", to compile the code.
Then: "./a.out", to execute it.
Output files are .dat.

-lm: calls the math library of your machine
-fopenmp: used for multithreading (enables pragma)

This codes can be adapted for different chemical models to generate random numbers for the initial conditions.
For this example, an Euler algorithm is shown. A 4th-order Runge-Kutta can easily implemented to replace Euler's scheme.

*/
/***********************************************************************************************/
//These are the libraries
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
// These are the files
FILE *fichierFront;
char filenameFront[200];
//Some definitions
double A,R,S,resA,resR,resS;
double ks,ksr,ka,kar,km,kn,knr;
double SumA,SumR,SumS;
int N,T,k,t,p,i, j,iprint;                  
double step;                        /*step = D / dx²*/
int tp = 0;
   double fp;
/*----------------------------------------------< End: declaration >---------------------------------------------------------------------------------------------*/
/*Main*/
int main()
{
  double L = 250.;    //size of the system                                 
  double dx = 0.5;  //spatial increment
  double Ttot = 200000.; //total simulation time
  double dt = 0.1; //time increment
  //Front position
  double fp,old_fp,fp_1,fp_2;
  old_fp = dx*((N-2)/2);
    //Number of integration steps and vector size
    N = (int)(L/dx);
    T = (int)(Ttot/dt);
    /*Concentration vectors */
    double CA[N]; 
    double CR[N];
    double CS[N];
    double CAnew[N];
    double CRnew[N];
    double CSnew[N];
    //Files
    sprintf(filenameFront,"Front_kdpi.dat");
    fichierFront = fopen(filenameFront, "w");
    //values for CDA parameter gamma
    double gama_sample[10] = {-0.5,-0.25,-0.1,-0.05,-0.01,0.01,0.05,0.1,0.25,0.5};
    //loop over gamma values
    int q;
    for ( q = 0; q < 10; q++)
    {
    //parameters and concentrations
    double gama = gama_sample[q];
    double pente = 0.;
    double eps= 1000.;
    double fac= 0.001;
    double D = 0.001;
    ks = fac*.5;
    ksr = ks/eps;
    ka = fac*1.;
    kar = ka/eps;
    km = fac*1.;
    kn = fac*0.1;
    knr = fac*0.5;
    A = 1.;
    R = 0.;
    S = 0.;
    //step
    step = 1/(dx*dx);
    /***********************************************************************************************/
    /*Initial conditions*/
    //S Domain
    for (i = 1 ; i <= (N-2)/2 ; i++)
    {
        CA[i] = A;
        CR[i] = (CA[i]*ka*kar - CA[i]*ka*km - kar*ksr + km*ksr - sqrt((kar - km)*(4*CA[i]*kar*kar*ks + (kar - km)*(-CA[i]*ka + ksr)*(-CA[i]*ka + ksr))))/(2*kar*(kar - km));
        CS[i] = (CA[i]*ka*kar - CA[i]*ka*km - kar*ksr + km*ksr + sqrt((kar - km)*(4*CA[i]*kar*kar*ks + (kar - km)*(-CA[i]*ka + ksr)*(-CA[i]*ka + ksr))))/(2*kar*(kar - km));
    }
    //R Domain
    for (i = 1 + (N-2)/2 ; i < N-1 ; i++)
    {
        CA[i] = A;
        CR[i] = (CA[i]*ka*kar - CA[i]*ka*km - kar*ksr + km*ksr + sqrt((kar - km)*(4*CA[i]*kar*kar*ks + (kar - km)*(-CA[i]*ka + ksr)*(-CA[i]*ka + ksr))))/(2*kar*(kar - km));
        CS[i] = (CA[i]*ka*kar - CA[i]*ka*km - kar*ksr + km*ksr - sqrt((kar - km)*(4*CA[i]*kar*kar*ks + (kar - km)*(-CA[i]*ka + ksr)*(-CA[i]*ka + ksr))))/(2*kar*(kar - km));
    }
    //Domains boundary
    int mid = (N-2)/2;
    CR[mid]=sqrt((CR[mid-1]-CR[mid+1])*(CR[mid-1]-CR[mid+1]));
    CS[mid]=sqrt((CS[mid-1]-CS[mid+1])*(CS[mid-1]-CS[mid+1]));
    /***********************************************************************************************/
    /*First integration loops --> Let the system relax from its initial condition*/
    k = 0;
    int Trelax = (int)(T/20.);
    for (t = 0; t <= Trelax; t++)
    {
      /*Boundary conditions*/
      CR[0] = (4.*CR[1] - CR[2])/3.;
      CS[0] = (4.*CS[1] - CS[2])/3.;;
      CR[N] = (4.*CR[N-1] - CR[N-2])/3.;
      CS[N] = (4.*CS[N-1] - CS[N-2])/3.;
        for (i = 1; i < N-1; i++)
        {
            resR = CR[i] + dt*(ks*CA[i]-ksr*CR[i]+ka*CA[i]*CR[i]-kar*CR[i]*CR[i]-km*CR[i]*CS[i]);
            resS = CS[i] + dt*(ks*CA[i]-ksr*CS[i]+ka*CA[i]*CS[i]-kar*CS[i]*CS[i]-km*CR[i]*CS[i]);
            CRnew[i] = D*dt*step*(1. + gama)*(CR[i+1] + CR[i-1] - 2.*CR[i]) + resR;
            CSnew[i] = D*dt*step*(1. - gama)*(CS[i+1] + CS[i-1] - 2.*CS[i]) + resS;
        }
        /* Update des concentrations */
        for (i = 1; i < N-1; i++)
        {
            CR[i] = CRnew[i];
            CS[i] = CSnew[i]; 
        }
    }

    /***********************************************************************************************/
    /*Second integration loops --> Track the velocity*/
    
    for(t = 0 ; t <= T ; t++)// Time loop 
    {
        /*Boundary conditions*/               
        CR[0] = (4.*CR[1] - CR[2])/3.;
        CS[0] = (4.*CS[1] - CS[2])/3.;;
        CR[N] = (4.*CR[N-1] - CR[N-2])/3.;
        CS[N] = (4.*CS[N-1] - CS[N-2])/3.;
        /*Initialization*/
        SumA = 0.;
        SumR = 0.;
        SumS = 0.;
        iprint = 1;
        for (i = 1; i < N-1; i++) // 1D spatial loop        
        {
            resR = CR[i] + dt*(ks*CA[i]-ksr*CR[i]+ka*CA[i]*CR[i]-kar*CR[i]*CR[i]-km*CR[i]*CS[i]);
            resS = CS[i] + dt*(ks*CA[i]-ksr*CS[i]+ka*CA[i]*CS[i]-kar*CS[i]*CS[i]-km*CR[i]*CS[i]);
            CRnew[i] = D*dt*step*(1. + gama)*(CR[i+1] + CR[i-1] - 2.*CR[i]) + resR;
            CSnew[i] = D*dt*step*(1. - gama)*(CS[i+1] + CS[i-1] - 2.*CS[i]) + resS;
            // Total concentrations
            SumR += CR[i];
            SumS += CS[i];
        }//End of spatial loops 

        /***********************************************************************************************/

        /* Front position and results printing*/
          int f = 1;
          while (CR[f] > CS[f]) //scan the system to find the front position
          {
            fp = f*dx;
            f += 1;
          } 
        if (t==T/4)
        {
          fp_1 = fp;
        }
        if (t==T)
        {
          fp_2 = fp;
        }
        old_fp = fp;
        /*Concentrations update*/
        for (i = 1; i < N-1; i++)
        {
            CR[i] = CRnew[i];
            CS[i] = CSnew[i]; 
        }
    }  //End of time loop
  // Printing of the front position in the result file
   fprintf(fichierFront,"%lf \t %lf \t %lf \n",gama,(fp_2-fp_1),(fp_2-fp_1)/(0.75*T*dt));
  }
    return 0;
}

/*------------------------------------------------------------------< Fin >---------------------------------------------------------------------------------------*/