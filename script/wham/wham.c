/*
   (C) Emilio Gallicchio, Rutgers University, 2005

   A program for WHAM (Weighted Histogram Analysis Method) analysis of 
   histograms with error analysis. Use of this program in published work
   should be acknowledged by citing the following paper:

   Gallicchio, E., M. Andrec, A.K. Felts, and R.M. Levy.  
   "T-WHAM, Replica Exchange, and Transition Paths." 
    J. Phys. Chem. B, (2005) in press. 

 */

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "ranlib.h"

#define MAXSTRING 128

/*extern int dirichlet_gen2bins(int nbins, double *p, double *counts, 
			      int bin0, int bin1);*/


int wham_error_analysis(int nbins, double *counts, double *p0,
			int nsim, double *simcounts, double **cu,
			double *avp, double *dp, int nsamples,
			int npasses, int no_zerocb, int verbose,
			FILE *samples_file){
  /* 
  Samples WHAM distribution and returns standard deviation of each
  population.

  P({p}|{n}) = [Prod_i (fi)^Ni] [Prod_j (pj)^nj]

    fi = Sum_j cij pj

  where:
  i: simulation index
  Ni: total count for simulation i
  j: bin index
  pj: probability of bin j
  nj: total count in bin j (from any simulation)
     
  Samples by rejection using proposals from Dirichlet distribution 
  */

  double *p = (double *)calloc(nbins,sizeof(double));
  double *av2p = (double *)calloc(nbins,sizeof(double));
  double *f = (double *)calloc(nsim,sizeof(double));
  double *fnew = (double *)calloc(nsim,sizeof(double));
  int i,j, jj, m;
  long moves = 0;
  long accepted = 0;
  int pass; 
  long seed1 = (long)time(NULL);
  long seed2 = (long)sqrt((double)seed1);
  int bin0, bin1;
  double p1_old, p0_old;
  double t;
  float csi;
  int sample;

  /* intialize starting probabilities and accumulators */
  for(j=0;j<nbins;j++){
    p[j] = p0[j];
    dp[j] = 0.0;
    avp[j] = 0.0;
  }

  /* calculate initial f's */
  for(i=0;i<nsim;i++){
    f[i] = 0.0;
    for(j=0;j<nbins;j++){
      f[i] += cu[i][j]*p0[j];
    }
  }

  /* init random number generator */
  setall(seed1,seed2);

  if(verbose) fprintf(stderr,"Starting sampling ...\n");

  for(sample=0;sample<nsamples;sample++){

    if(verbose) fprintf(stderr,"Begin sample %d\n",sample+1);

    for(pass=0;pass<npasses;pass++){

      if(verbose) fprintf(stderr,"Begin pass %d\n",pass+1);

      for(bin1=0;bin1<nbins;bin1++){
	
	if(no_zerocb && counts[bin1] <= 0) continue;
	
	moves += 1;
	
	/* select another bin at random */
	bin0 = bin1;
	while(bin0 == bin1 || (no_zerocb && counts[bin0] <= 0)){
	  bin0 = (int)ignuin(0,nbins-1);
	}
	p0_old = p[bin0];
	p1_old = p[bin1];
	/*if(dirichlet_gen2bins(nbins,p,counts,bin0,bin1)<0){
	  fprintf(stderr,"wham_error_analysis: error in dirichlet_gen2bins().\n");
	  return -1;
	}*/
	
	/* recompute f's */
	for(i=0;i<nsim;i++){
	  fnew[i] = f[i] + 
	    cu[i][bin0]*(p[bin0]-p0_old) + cu[i][bin1]*(p[bin1]-p1_old);
	}
	t = 1.0;
	for(i=0;i<nsim;i++){
	  t *= pow(f[i]/fnew[i],simcounts[i]);
	}
	csi =  genunf(0.,1.); /* random numb  0 <= x <= 1 */
	
	if(t>(double)csi){
	  /* accept */
	  accepted += 1;
	  for(i=0;i<nsim;i++){
	    f[i] = fnew[i];
	  }
	}else{
	  /* reject */
	  p[bin0] = p0_old;
	  p[bin1] = p1_old;
	}
      }

      if(verbose) fprintf(stderr,"End of pass %d\n",pass+1);
    }

    for(j=0;j<nbins;j++){
      t = p[j];
      avp[j] += t;
      av2p[j] += t*t;
      if(samples_file){
	fprintf(samples_file,"%16.8e\n",t);
      }
    }
    if(samples_file) fflush(samples_file);

    if(verbose) fprintf(stderr,"End of sample %d\n",sample+1);

  }
  

  /* print acceptance ratio */
  fprintf(stderr,"Error analysis acceptance ratio: %10.3f\n",
	  (double)accepted/(double)moves);

  /* calculate standard deviations */
  for(j=0;j<nbins;j++){
    avp[j] /= (double)nsamples;
    av2p[j] /= (double)nsamples;
    dp[j] = sqrt((av2p[j]-avp[j]*avp[j]));
  }

  free(p);
  free(av2p);
  free(f);
  free(fnew);
  return 1;
}

void print_usage(void){
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"wham [-h] [-v] -c <counts_file> -s <simcount_file>\n-u <umbrella_potential_file>\n [-r <rmsd_convergence_value>] [-x <max_iterations>] \n[-e [-o <error_analysis_file>] [-n <nsamples>] [-p <passes_per_sample>]\n[-z] ]\n");
}


int main(int argc, char *argv[]){

  /* usage:
     wham -c <counts_file> -s <simcount_file> -u <umbrella_potential_file> 
      [-r <rmsd_convergence_value>] [-x <max_iterations>] 
      [-e [-o <error_analysis_file>] [-n <nsamples>] [-p <passes_per_sample>] 
          [-z] ]
      
      if -e is given do error analysis
      if -o is given dump <nsamples> alternative probability distributions
        in file <error_analysis_file>.
      if -z is given do zero counts bin sampling.

      Defaults:
       <rmsd_convergence_value> = 0.1
       <max_iterations> = 1000
       <nsamples> = 128
       <passes_per_sample> = 10
  */

  /*
     p(xj) = N(xj)/sum_i ni exp(fi) cij

     exp(-fi) = sum_j cij p(xj)

     where:
       p(xj): unbiased distribution at xj
       N(xj): total count at xj from all simulations (from counts_file)
       ni: total samples from simulation i (from totalcount_file)
       cij: umbrella potential in bin j of simulation i 
            (from umbrella_potential_file)

	    cij = exp[-(betai-beta)Ej] exp(-betai wi(x_j))

	    where betai is temperature of simulation of i
	    and wi(x_j) is umbrella potential at bin j from simulation i.
	    If the temperatue is different in each simulation it is assumed
	    that the potential energy E is one of the binned quantities. 
	    Therefore Ej is potential energy in bin j. If the temperature
	    is the same in all simulations the first term in the r.h.s of
	    the equation for cij is 1 (betai = beta).

  */

  int nsim; /* number of simulations */
  int nbins; /* number of bins */
  char counts_file_name[MAXSTRING] = "";
  char simcounts_file_name[MAXSTRING] = "";
  char umbrella_file_name[MAXSTRING] = "";
  char samples_file_name[MAXSTRING] = "";
  FILE *counts_file, *simcounts_file, *umbrella_file, *samples_file;
  char *ch;
  char line[MAXSTRING];
  double *counts;
  double *simcounts;
  double **cu;
  int i,j;
  double *u;
  double *f;
  double *p, *dp, *avp;
  double *t;
  double rms;
  double rmscut = 0.1;
  int maxiter = 1000;
  int nsamples = 128;
  int npasses = 10;
  int verbose = 0;
  int iter, l;
  double norm;
  int error_analysis = 0;
  int no_zerocb = 1;
  extern char *optarg;
  extern int optind, opterr, optopt;
  int copt;

  /* parse options */
  while((copt = getopt(argc, argv, "c:s:u:r:x:eo:n:p:zhv")) != -1){
    switch (copt){

    case '?':
      fprintf(stderr,"Warning: Unrecognized option\n");
      print_usage();
      break;
    case 'h':
      print_usage();
      exit(1);
      break;
    case ':':
      fprintf(stderr,"Warning: option requires a parameter\n");
      print_usage();
      break;
    case 'c':
      strcpy(counts_file_name,optarg);
      break;
    case 's':
      strcpy(simcounts_file_name,optarg);
      break;
    case 'u':
      strcpy(umbrella_file_name,optarg);
      break;
    case 'r':
      if(sscanf(optarg,"%lf",&rmscut) != 1){
	fprintf(stderr,"Warning: error reading rmsd convergence value\n");
	print_usage();
      }
      break;
    case 'x':
      if(sscanf(optarg,"%d",&maxiter) != 1){
	fprintf(stderr,"Warning: error reading max number of iterations\n");
	print_usage();
      }
      break;
    case 'e':
      error_analysis = 1;
      break;
    case 'o':
      strcpy(samples_file_name,optarg);
      break;
    case 'n':
      if(sscanf(optarg,"%d",&nsamples) != 1){
	fprintf(stderr,"Warning: error reading number of samples\n");
	print_usage();
      }
      break;
    case 'p':
      if(sscanf(optarg,"%d",&npasses) != 1){
	fprintf(stderr,"Warning: error reading number of passes per sample\n");
	print_usage();
      }
      break;
    case 'z':
      no_zerocb = 0;
      break;
    case 'v':
      verbose = 1;
      break;
    default:
      break;
    }
  }

  if(strlen(counts_file_name) == 0){
    fprintf(stderr,"Error: counts file is not defined\n");
    print_usage();
    return 0;
  }
  if(strlen(simcounts_file_name) == 0){
    fprintf(stderr,"Error: simulation counts file is not defined\n");
    print_usage();
    return 0;
  }
  if(strlen(umbrella_file_name) == 0){
    fprintf(stderr,"Error: umbrella file is not defined\n");
    print_usage();
    return 0;
  }

  /* read counts file */
  if((counts_file = fopen(counts_file_name, "r")) == NULL){
    fprintf(stderr,"wham: cannot open file %s for reading.\n",
            counts_file_name);
    return(-1);
  }
  nbins = 0;
  while((ch = (char *)fgets(line, MAXSTRING, counts_file)) != (char *)NULL){
    nbins += 1;
  }
  counts = (double *)calloc(nbins,sizeof(double));
  if(!counts){
    fprintf(stderr,"wham: cannot allocate counts array (%d doubles)\n",
	    nbins);
    return(-1);
  }  
  rewind(counts_file);
  for(j=0;j<nbins;j++){
    if(fscanf(counts_file, "%lf", &(counts[j])) != 1){
      fprintf(stderr,"wham: problem reading file %s.\n",
	      counts_file_name);
      return(-1);
    }
  }
  fclose(counts_file);

  /* read simulation counts file */
  if((simcounts_file = fopen(simcounts_file_name, "r")) == NULL){
    fprintf(stderr,"wham: cannot open file %s for reading.\n",
            simcounts_file_name);
    return(-1);
  }
  nsim = 0;
  while((ch = (char *)fgets(line, MAXSTRING, simcounts_file)) != (char *)NULL){
    nsim += 1;
  }
  simcounts = (double *)calloc(nsim,sizeof(double));
  if(!simcounts){
    fprintf(stderr,"wham: cannot allocate simcounts array (%d doubles)\n",
	    nsim);
    return(-1);
  }  
  rewind(simcounts_file);
  for(i=0;i<nsim;i++){
    if(fscanf(simcounts_file, "%lf", &(simcounts[i])) != 1){
      fprintf(stderr,"wham: problem reading file %s.\n",
	      simcounts_file_name);
      return(-1);
    }
  }
  fclose(simcounts_file);

  /* read umbrella potential file. Matrix in row order. Each row corresponds
     to a simulation and columns are bins. */
  cu = (double **)calloc(nsim,sizeof(double *));
  if(!cu){
    fprintf(stderr,"wham: cannot allocate cu array (%d double *)\n",
	    nsim);
    return(-1);
  }
  for(i=0;i<nsim;i++){
    cu[i] = (double *)calloc(nbins,sizeof(double));
    if(!cu[i]){
      fprintf(stderr,"wham: cannot allocate cu[%d] array (%d double *)\n",
	      i,nbins);
      return(-1);
    }
  }
  if((umbrella_file = fopen(umbrella_file_name, "r")) == NULL){
    fprintf(stderr,"wham: cannot open file %s for reading.\n",
            umbrella_file_name);
    return(-1);
  }
  l = 0;
  for(i=0;i<nsim;i++){
    for(j=0;j<nbins;j++){
      l += 1;
      if(fscanf(umbrella_file, "%lf", &(cu[i][j])) != 1){
	fprintf(stderr,"wham: problem reading file %s at line %d.\n",
		umbrella_file_name, l);
	return(-1);
      }
    }
  }
  fclose(umbrella_file);

  /* start setting fi = 1 */
  f = (double *)calloc(nsim,sizeof(double));
  u = (double *)calloc(nsim,sizeof(double));
  p = (double *)calloc(nbins,sizeof(double));
  t = (double *)calloc(nbins,sizeof(double));
  if(!f || !u || !p || !t){
    fprintf(stderr,"wham: cannot allocate work arrays (%d doubles)\n", 
	    2.*nsim+2*nbins);
    return(-1);
  }
  for(i=0;i<nsim;i++){
    f[i] = 1.0;
    u[i] = simcounts[i]/f[i];
  }

  /* calculate p until rms is satisfied */
  rms = 2.*rmscut;
  iter = 1;
  while(rms >= rmscut && iter < maxiter){

    for(j=0;j<nbins;j++){
      t[j] = 0.0;
      for(i=0;i<nsim;i++){
	t[j] += u[i]*cu[i][j];
      }
      p[j] = counts[j]/t[j];
    }

    /* normalize p */
    norm = 0.0;
    for(j=0;j<nbins;j++){
      norm += p[j];
    }
    for(j=0;j<nbins;j++){
      p[j] /= norm;
    }

    /* recalculates f's */
    for(i=0;i<nsim;i++){
      f[i] = 0.0;
      for(j=0;j<nbins;j++){
	f[i] += cu[i][j]*p[j];
      }
      u[i] = simcounts[i]/f[i];
    }

    /* rmsd from data */
    rms = 0.0;
    for(j=0;j<nbins;j++){
      t[j] = 0.0;
      for(i=0;i<nsim;i++){
        t[j] += (simcounts[i]/f[i])*cu[i][j]*p[j];
      }
      rms += pow(counts[j]-t[j],2);
    }
    rms = sqrt(rms/(double)nbins);

    fprintf(stderr,"Iteration: %d. RMSD= %lf\n", iter, rms);
    iter += 1;
  }

  if(rms<rmscut){
    fprintf(stderr,"Convergence reached.\n");
  }else{
    fprintf(stderr,"Maximum number of iterations exceeded.\n");
  }
  
  fprintf(stderr,"Simulation   -log(f)\n");
  for(i=0;i<nsim;i++){
    fprintf(stderr,"%d %16.8e %.8f\n",i,-log(f[i]), f[i]);
  }


  if(error_analysis){
    dp = (double *)calloc(nbins,sizeof(double));
    avp = (double *)calloc(nbins,sizeof(double));
    if(strlen(samples_file_name) == 0){
      samples_file = NULL;
    }else{
      samples_file = fopen(samples_file_name,"w");
      if(!samples_file){
	fprintf(stderr,"Unable to open file %d for writing\n",
		samples_file_name);
	return 0;
      }
    }
    wham_error_analysis(nbins,counts,p,nsim,simcounts,cu,avp,dp,
			nsamples,npasses,no_zerocb,verbose,samples_file);
  }

  if(error_analysis){
    for(j=0;j<nbins;j++){
      printf("%16.8e %16.8e %16.8e\n",p[j],avp[j],dp[j]);
    }
  }else{
    for(j=0;j<nbins;j++){
      printf("%16.8e\n",p[j]);
    }
  }

}
