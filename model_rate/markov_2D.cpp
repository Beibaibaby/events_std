#include "mex.h"
#include <math.h>
#include <random>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>

// Constants for ran2 random number generator. From Numerical Recipes.
#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

// define some convenient functions
double ran2(long *idum) ;    // Ran2 uniform random deviate generator from Numerical Recipes
double sigtransfer(double x); //sigmoidal transfer function
double twoscale_pf(double rate, double tr, double td, double thresh, double mag, double beta);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
        
    //simulation inputs
    double *times, *params, *ee_params, *ie_params;
    double *N, *n0, *rng_stuff;
    
    //sim outputs
    double *return_t, *return_ve, *return_vi, *return_p, *return_p2;
    
    //working arrays
    double *ve, *vi, *plastics, *intensity, *probs, *ctrarray, *p2;

    //get parameter values into objects
    times      = mxGetPr(prhs[0]);
    N          = mxGetPr(prhs[1]);
    n0         = mxGetPr(prhs[2]);
    params     = mxGetPr(prhs[3]);
    ee_params  = mxGetPr(prhs[4]);
    ie_params  = mxGetPr(prhs[5]);
    rng_stuff  = mxGetPr(prhs[6]);
    
    
    //define parameter names
    double jee, jie, jei, jii, Ie, Ii, taui;
    double eps, gamma;
    double trie, tdie, threshie, magie, slie;
    double tree, tdee, threshee, magee, slee;
    double tlim, ttoss,dt;
    int ne, ni;
    double pee, pie, ree, rie, rei, rii, cee, cie, alphee, alphie;
    
    //various counters
    int point,s,j;
    s=0;
    j=0;
    double curtime, lasttime;
    int curpos, lastpos;

    //various working variables
    double smallcheck;
    double Ne, Ni;
    double ex_cur, in_cur, deltat, choice, pnew;

    //random number seed
    long idum;
    double *seed;
    
    //assign parameter values
    jee  = params[0];
    jei  = params[1];
    jie  = params[2];
    jii  = params[3];
    Ie   = params[4];
    Ii   = params[5];
    taui = params[6];
    
    tlim = times[0];
    ttoss = times[1];
    dt = times[2];
    
    ne = (int) n0[0];
    ni = (int) n0[1];

    tree     = ee_params[0];
    tdee     = ee_params[1];
    threshee = ee_params[2];
    magee    = ee_params[3];
    slee     = ee_params[4];

    trie     = ie_params[0];
    tdie     = ie_params[1];
    threshie = ie_params[2];
    magie    = ie_params[3];
    slie     = ie_params[4];

    //random number stuff
    /* somewhere in Matlab before you call this function, set seed = -round(10000*rand)
        and pass it into one of the parameter arrays
    */
    seed = &rng_stuff[0]; 
    idum=seed[0]<0?(long)seed[0]:-(long)seed[0]; // for random number generator

    //allocate memory for working arrays
    int maxsize;
    maxsize = ceil((tlim-ttoss)/dt);

    ve = (double *) mxMalloc(maxsize * sizeof(double));         //excitatory population
    vi = (double *) mxMalloc(maxsize * sizeof(double));         //inhibitory population
    plastics = (double *) mxMalloc(maxsize * sizeof(double));   //pee values
    p2 = (double *) mxMalloc(maxsize * sizeof(double));         //pie values
    ctrarray = (double *) mxMalloc(maxsize * sizeof(double));   //for discretizing into temporally uniform bins
    intensity = (double *) mxMalloc(4 * sizeof(double));        //for gillespie
    probs = (double *) mxMalloc(4 * sizeof(double));            //for gillespie
    //initialize the arrays to zeros
    for (j=0;j<4;j++){
        intensity[j]=0;
        probs[j]=0;
    }
    for (j=0;j<maxsize;j++){
        ve[j]=0;
        vi[j]=0;
        plastics[j]=0;
        p2[j]=0;
        ctrarray[j]=0;
    }

    //initialize various and sundry values
    smallcheck = pow(10,-15);
    Ne = N[0];
    Ni = N[1];
    ve[0] = n0[0]/Ne;
    vi[0] = n0[1]/Ni;

    //initialize plasticity onto the initial condition
    pee = twoscale_pf(ne/Ne,tree,tdee,threshee,magee,slee);
    alphee = 1/(tree*pee);
    pie = twoscale_pf(ne/Ne,trie,tdie,threshie,magie,slie);
    alphie = 1/(trie*pie);

    plastics[0]=pee;
    p2[0]=pie;
    ctrarray[0]=1;

    //initialize rates
    ree = (ne/Ne)*pee;
    rei = ni/Ni;
    rie = (ne/Ne)*pie;
    rii = ni/Ni;

    //zero out starting positions for the loop
    curtime=0;
    lasttime=0;
    curpos=0;
    lastpos=0;


    while (curtime < tlim){
        //set current position based on current time
        curpos = floor((curtime-ttoss)/dt);

        //if ttoss>0, check that we aren't in the initial transient to ignore
        if (curtime<ttoss) curpos=0;

        //initialize total intensity to zero
        s=0;

        //find synapse relaxation constant        
        cee = twoscale_pf(ne/Ne,tree,tdee,threshee,magee,slee);
        alphee = 1/(tree*cee);

        cie = twoscale_pf(ne/Ne,trie,tdie,threshie,magie,slie);
        alphie = 1/(trie*cie);
        
        //set input currents
        ex_cur = jee*ree - jei*rei + Ie;
        in_cur = jie*rie - jii*rii + Ii;
        
        //initialize intensities for gillespie
        intensity[0] = (double) ne;
        intensity[1] = (double) Ne*sigtransfer(sqrt(Ne)*ex_cur);
        intensity[2] = (double) (1/taui)*ni;
        intensity[3] = (double) (1/taui)*Ni*sigtransfer(sqrt(Ni)*in_cur);

        //truncate to ensure we can't get bigger than N
        if (ne >= Ne){
            intensity[1]=0;
        }
        if (ni >= Ni){
            intensity[3]=0;
        }

        //find total intensity for normalization        
        for (j=0;j<4;j++){
            s += intensity[j];
        }
            
        //when low intensity, ensure no division by zero
        if (s<smallcheck){
            printf("s is very small\n");
            s+=pow(10,-6);
        }

        //set probability intervals        
        probs[0] = intensity[0]/s;
        for (j=1;j<4;j++){
            probs[j] = probs[j-1]+(intensity[j]/s);
        }
        
        //find time of next event
        deltat = -log(ran2(&idum))/s;
        
        //let synapse relax to new value
        pnew = cee-(cee-pee)*(exp(-alphee*deltat));
        pee = pnew;
               
        pnew = cie-(cie-pie)*(exp(-alphie*deltat));
        pie = pnew;

        //choose which event happens        
        choice=ran2(&idum);
        if (choice < probs[0]){
            ne = ne-1;
        }else if (choice < probs[1]){
            ne = ne+1;
        }else if (choice < probs[2]){
            ni = ni-1;
        }else{
            ni = ni+1;
        }
        
        //update postsynaptic instantaneous rate
        ree = pee*ne/Ne;
        rie = pie*ne/Ne;
        rei = ni/Ni;
        rii = ni/Ni;
        
        //update the output vectors
        ve[curpos] = ve[curpos] + ne/Ne;
        vi[curpos] = vi[curpos] + ni/Ni;
        plastics[curpos] = plastics[curpos]+pee;
        p2[curpos] = p2[curpos]+pie;
        ctrarray[curpos] = ctrarray[curpos]+1;

        //make sure we haven't skipped over a bin in the discretization
        if ((curpos - lastpos) > 1){
            //if we skipped a cell, fill in all the intervening cells as though
            //everything was constant
            for (j=lastpos+1;j<curpos;j++){
                ve[j] = ve[lastpos];
                vi[j] = vi[lastpos];
                plastics[j] = plastics[lastpos];
                p2[j] = p2[lastpos];
                ctrarray[j] = ctrarray[lastpos];
            }
        }

        curtime = curtime + deltat;
        lastpos = curpos;
    }
    

    //now we assign the return variables for matlab to deal with
    plhs[0] = mxCreateDoubleMatrix( curpos, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix( curpos, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix( curpos, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix( curpos, 1, mxREAL);
    plhs[4] = mxCreateDoubleMatrix( curpos, 1, mxREAL);

    return_t  = mxGetPr(plhs[0]);
    return_ve = mxGetPr(plhs[1]);
    return_vi = mxGetPr(plhs[2]);
    return_p  = mxGetPr(plhs[3]);
    return_p2 = mxGetPr(plhs[4]);

    for (j=0;j<curpos;j++){
        return_t[j]  = j*dt;
        return_ve[j] = ve[j]/ctrarray[j];       //normalize for the averaging
        return_vi[j] = vi[j]/ctrarray[j];
        return_p[j]  = plastics[j]/ctrarray[j];
        return_p2[j] = p2[j]/ctrarray[j];
    }

    //free up memory
    mxFree(ve);
    mxFree(vi);
    mxFree(plastics);
    mxFree(p2);
    mxFree(ctrarray);
    mxFree(intensity);
    mxFree(probs);

    return;
}
        
//helpful functions
//Numerical Recipes ran2
double ran2(long *idum){
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    float temp;

    if (*idum <= 0) {
        if (-(*idum) < 1) 
            *idum=1;
        else 
            *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0){
                *idum += IM1;}
            if (j < NTAB){
                iv[j] = *idum;}
        }
        iy=iv[0];
    }

    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
}


//sigmoidal transfer function
double sigtransfer(double x){
    double f;
    f = 1.0/(1.0+exp(-x));
    return f;
}

//for plasticity
double twoscale_pf(double rate, double tr, double td, double thresh, double mag, double beta){
    double a, pf;
    a =  (double) mag/(1+exp(-beta*(rate-thresh)));
    pf = (double) td/(td + tr*a);
    return pf;
}

