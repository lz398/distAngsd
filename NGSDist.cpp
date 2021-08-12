/*
 *
 * ngsDist - NGS data individual inbreeding coefficients estimation.
 * Copyright (C) 2012  Filipe G. Vieira
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
*/

#include "ngsDist.hpp"
#include "emOptim2.cpp"

char const* version = "1.0.9";

void rnd_map_data(params *pars, uint64_t n_blocks);

int main (int argc, char** argv) {
  /////////////////////
  // Parse Arguments //
  /////////////////////
  uint64_t n_lines, n_fields;
  params* pars = new params;
  init_pars(pars);
  parse_cmd_args(pars, argc, argv);



  ///////////////////////
  // Adjust Parameters //
  ///////////////////////
  // Calculate total number of combinations
  uint64_t n_comb = (pow(pars->n_ind, 2) - pars->n_ind) / 2;
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Analysis will be run in %lu combinations\n", n_comb);
  // Adjust thread number to combinations
  if(n_comb < pars->n_threads){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Fewer combinations (%ld) than threads (%d). Reducing the number of threads...\n", n_comb, pars->n_threads);
    pars->n_threads = n_comb;
  }

  // If input are genotypes (either called ot not) assume independence between genotypes (faster)
  if(!pars->in_probs && !pars->indep_geno){
    fprintf(stderr, "==> Using faster algorithm (assuming independence of genotypes) since input are genotypes!\n");
    pars->indep_geno = true;
  }
  else if(pars->call_geno && !pars->indep_geno){
    fprintf(stderr, "==> Using faster algorithm (assuming independence of genotypes) since calling genotypes!\n");
    pars->indep_geno = true;
  }
  else if(pars->indep_geno)
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Using faster algorithm (assuming independence of genotypes)!\n");



  ///////////////////////
  // Check input files //
  ///////////////////////
  // Get file total size
  if( strcmp(pars->in_geno, "-") == 0 ){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Reading from STDIN (BINARY)\n");
    pars->in_bin = true;
  }else{
    struct stat st;
    if( stat(pars->in_geno, &st) != 0 )
      error(__FUNCTION__, "cannot check GENO file size!");

    if( strcmp(strrchr(pars->in_geno, '.'), ".gz") == 0 ){
      if(pars->verbose >= 1)
    fprintf(stderr, "==> GZIP input file (never BINARY)\n");
      pars->in_bin = false;
    }else{
      if(pars->verbose >= 1)
    fprintf(stderr, "==> BINARY input file\n");
      pars->in_bin = true;
      pars->in_probs = true;

      if( pars->n_sites != st.st_size/sizeof(double)/pars->n_ind/N_GENO )
    error(__FUNCTION__, "invalid/corrupt genotype input file!");
    }
  }
  


  ////////////////////////////
  // Prepare initial values //
  ////////////////////////////
  // Read labels file
  if(pars->in_labels){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Reading labels\n");
    n_lines = read_file(pars->in_labels, &pars->ind_labels, (pars->in_labels_header ? 1 : 0));
    if(n_lines != pars->n_ind)
      error(__FUNCTION__, "invalid LABELS file!");

    // Fix labels...
    char* ptr;
    for(uint64_t i = 0; i < pars->n_ind; i++){
      ptr = strchr(pars->ind_labels[i], '\t');
      if(ptr != NULL)
    *ptr = '\0';
    }
  }else{
    pars->ind_labels = init_ptr(pars->n_ind, BUFF_LEN, (const char*) "Ind_#");
    // Tweak initiation value (replace # by number)
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* pch = strchr(pars->ind_labels[i], '#');
      if(pch)
    sprintf(pch, "%lu", i);
    }
  }
  if(pars->verbose >= 4)
    for(uint64_t i = 0; i < pars->n_ind; i++)
      fprintf(stderr, "%s\n", pars->ind_labels[i]);



  // Read POS file
  if(pars->in_pos){
    if(pars->verbose >= 1)
      fprintf(stderr, "==> Reading positions file\n");

    // Allocate memory
    pars->pos = init_ptr(pars->n_sites, 0, 0, (const char*) '\0');

    // Read and split file
    read_split(pars->in_pos, pars->pos, &n_lines, &n_fields, (pars->in_pos_header ? 1 : 0));

    if(n_lines != pars->n_sites || n_fields < 2)
      error(__FUNCTION__, "invalid POS file!");

    if(pars->verbose >= 4)
      for(uint64_t s = 0; s < pars->n_sites; s++)
    fprintf(stderr, "%s\t%s\t%s\t%s\n", pars->pos[s][0], pars->pos[s][1], pars->pos[s][2], pars->pos[s][3]);
  }



  // Read GENO file
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Reading genotype data\n");
  pars->in_geno_lkl = read_geno(pars->in_geno, pars->in_bin, pars->in_probs, &pars->in_logscale, pars->n_ind, pars->n_sites);



  // Make copy of GLs in case of bootstrap
  pars->geno_lkl = init_ptr(pars->n_ind, pars->n_sites, 0, -INF);
  for(uint64_t i = 0; i < pars->n_ind; i++)
    memcpy(pars->geno_lkl[i], pars->in_geno_lkl[i], (pars->n_sites)*sizeof(double*));

  for(uint64_t i = 0; i < pars->n_ind; i++)
    for(uint64_t s = 0; s < pars->n_sites; s++){
      // Call genotypes
      if(pars->call_geno)
    call_geno(pars->in_geno_lkl[i][s], N_GENO, pars->in_logscale, pars->N_thresh, pars->call_thresh, 0);

      // Convert space - from now on, all in NORMAL space!
      if(pars->in_logscale)
    conv_space(pars->in_geno_lkl[i][s], N_GENO, exp);
    }

  // Initialize random number generator
  if(pars->verbose >= 2)
    fprintf(stderr, "==> Setting seed for random number generator\n");
  pars->rnd_gen = gsl_rng_alloc(gsl_rng_taus);
  gsl_rng_set(pars->rnd_gen, pars->seed);

 

  //////////////////////
  // Open output file //
  //////////////////////
  FILE* out_fh = fopen(pars->out, "w");
  if(out_fh == NULL)
    error(__FUNCTION__, "cannot open output file!");
 


  /////////////////////
  // Prepare threads //
  /////////////////////
  // Create pthread structure array
  pth_struct* pth = new pth_struct[n_comb];
  // Initialize pars and distance matrix pointers
  uint64_t comb_id = 0;
  double** dist_matrix = init_ptr(pars->n_ind, pars->n_ind, 0.0);
  for(comb_id = 0; comb_id < n_comb; comb_id++){
    pth[comb_id].pars = pars;
    pth[comb_id].dist_matrix = dist_matrix;
  }

  // Create threadpool
  if( (pars->thread_pool = threadpool_create(pars->n_threads, n_comb, 0)) == NULL )
    error(__FUNCTION__, "failed to create thread pool!");



  //////////////////
  // Analyze Data //
  //////////////////
  fflush(stdout);
  // Loop for bootstrap analyses
  for(uint64_t rep = 0; rep <= pars->n_boot_rep; rep++){
    comb_id = 0;

    if(pars->verbose >= 1){
      if(rep == 0)
    // Full dataset analyses
    fprintf(stderr, "==> Analyzing full dataset...\n");
      else
    // Bootstrap analyses
    fprintf(stderr, "==> Bootstrap replicate # %lu ...\n", rep);
    }


    // Map from in_geno_lkl data to geno_lkl
    if(pars->verbose >= 2)
      fprintf(stderr, "> Mapping positions...\n");

    // If not rep 0, then adjust number of sites and random sample blocks from original dataset
    if(rep > 0){
      pars->n_sites -= pars->n_sites % pars->boot_block_size;
      rnd_map_data(pars, pars->n_sites/pars->boot_block_size);
    }

    // Calculate pairwise genetic distances
    if(pars->verbose >= 2)
      fprintf(stderr, "> Calculating pairwise genetic distances...\n");

    for(uint64_t i1 = 0; i1 < pars->n_ind; i1++)
      for(uint64_t i2 = i1+1; i2 < pars->n_ind; i2++){
    // Set which individuals to analyze
    pth[comb_id].i1 = i1;
    pth[comb_id].i2 = i2;

    // Add task to thread pool
    int ret = threadpool_add(pars->thread_pool, gen_dist_slave, (void*) &pth[comb_id++], 0);
    if(ret == -1)
          error(__FUNCTION__, "invalid thread pool!");
    else if(ret == -2)
          error(__FUNCTION__, "thread pool lock failure!");
    else if(ret == -3)
          error(__FUNCTION__, "queue full!");
    else if(ret == -4)
          error(__FUNCTION__, "thread pool is shutting down!");
    else if(ret == -5)
          error(__FUNCTION__, "thread failure!");
      }


  
    //////////////////////////
    // Wait for all threads //
    //////////////////////////
    threadpool_wait(pars->thread_pool, 0.1);

    if(n_comb != comb_id)
      error(__FUNCTION__, "some combinations are missing!");



    ///////////////////////////
    // Print Distance Matrix //
    ///////////////////////////
    if(pars->verbose >= 2)
      fprintf(stderr, "> Printing distance matrix\n");

    fprintf(out_fh, "\n%lu\n", pars->n_ind);
    for(uint64_t i = 0; i < pars->n_ind; i++){
      char* buf = join(dist_matrix[i], pars->n_ind, "\t");
      fprintf(out_fh, "%s\t%s\n", pars->ind_labels[i], buf);
      delete [] buf;
    }

  }

  threadpool_wait(pars->thread_pool, 0.1);
  if(threadpool_destroy(pars->thread_pool, threadpool_graceful) != 0)
    error(__FUNCTION__, "cannot free thread pool!");

  fclose(out_fh);



  /////////////////
  // Free Memory //
  /////////////////
  if(pars->verbose >= 1)
    fprintf(stderr, "==> Freeing memory...\n");

  free_ptr((void**) dist_matrix, pars->n_ind);
  delete [] pth;
  // pars struct
  free_ptr((void***) pars->in_geno_lkl, pars->n_ind, pars->n_sites);
  free_ptr((void**) pars->geno_lkl, pars->n_ind);
  free_ptr((void**) pars->ind_labels, pars->n_ind);
  //free_ptr((void*) pars->in_geno);
  gsl_rng_free(pars->rnd_gen);

  if(pars->verbose >= 1)
    fprintf(stderr, "Done!\n");

  delete pars;

  return 0;
}




double gen_dist(params *p, uint64_t i1, uint64_t i2){
  uint64_t cnt = 0;
  double dist = 0;

  Matrix<double> GL1 = alloc(1,N_GENO);
  Matrix<double> GL2 = alloc(1,N_GENO);
  int dim = GL1.y*GL2.y;

  for(uint64_t s = 0; s < p->n_sites; s++){
    // Skip missing data
    if( p->pairwise_del &&
    (miss_data(p->geno_lkl[i1][s]) ||
     miss_data(p->geno_lkl[i2][s])) )
      continue;

    double* sfs = init_ptr(dim, (double) 1/dim);
    GL1.mat[0][0] = p->geno_lkl[i1][s][0];
    GL1.mat[0][1] = p->geno_lkl[i1][s][1];
    GL1.mat[0][2] = p->geno_lkl[i1][s][2];
    GL2.mat[0][0] = p->geno_lkl[i2][s][0];
    GL2.mat[0][1] = p->geno_lkl[i2][s][1];
    GL2.mat[0][2] = p->geno_lkl[i2][s][2];

    if(!p->indep_geno)
      em2(sfs, &GL1, &GL2, 0.001, 50, dim);

    for(uint64_t g1 = 0; g1 < N_GENO; g1++)
      for(uint64_t g2 = 0; g2 < N_GENO; g2++){
    dist += p->score[g1][g2] * (p->indep_geno ? p->geno_lkl[i1][s][g1]*p->geno_lkl[i2][s][g2] : sfs[3*g1+g2]);

    if(p->verbose >= 9)
      fprintf(stderr, "%lu\t%lu <-> %lu\t%lu - %lu\t%f\t%f\n", s, i1, i2, g1, g2, p->geno_lkl[i1][s][g1]*p->geno_lkl[i2][s][g2], sfs[3*g1+g2]);
      }

    if(p->verbose >= 8)
      fprintf(stderr, "Cumulative distance between %s (ind %lu) and %s (ind %lu) at site %lu: %f\n", p->ind_labels[i1], i1, p->ind_labels[i2], i2, s, dist);

    cnt++;
    free_ptr(sfs);
  }

  if(p->verbose >=3)
    fprintf(stderr, "\tDistance of %f from %lu valid sites (%f) between %s (ind %lu) and %s (ind %lu)!\n", dist, cnt, dist/(double) cnt, p->ind_labels[i1], i1, p->ind_labels[i2], i2);

  dalloc(GL1, 1);
  dalloc(GL2, 1);

  if(p->tot_sites > 0)
    cnt = p->tot_sites;

  // Calculates raw distance
  dist /= (double) cnt;
  // Evolutionary model
  if(p->evol_model == 0) {
    // Logarithmic transformation (log(1-d)) to make distance additive (assuming constant Ne) and avoid violating minimum evolution and NJ assumptions.
    dist = -log(1-dist);
  } else if(p->evol_model == 1) {
    // JC69
    dist = -log(1 - (dist * 4/3)) * 3/4;
  } else if(p->evol_model == 2) {
    // K80
    error(__FUNCTION__, "K80 model not yet supported");
  } else if(p->evol_model == 3) {
    // F81
    error(__FUNCTION__, "F81 model not yet supported");
  } else if(p->evol_model == 4) {
    // HKY85
    error(__FUNCTION__, "HKY85 model not yet supported");
  } else if(p->evol_model == 5) {
    // TN93
    error(__FUNCTION__, "TN93 model not yet supported");
  } else {
    error(__FUNCTION__, "invalid evolutionary model specified!");
  }

  return dist;
}



void gen_dist_slave(void *pth){
  pth_struct* p = (pth_struct*) pth;

  p->dist_matrix[p->i1][p->i2] = p->dist_matrix[p->i2][p->i1] = gen_dist(p->pars, p->i1, p->i2);
}



void rnd_map_data(params *pars, uint64_t n_blocks){
  uint64_t block, rnd_block;
  uint64_t block_s, rnd_block_s;

  // For each of the original blocks
  for(block = 0; block < n_blocks; block++){
    // Pick a random one to replace it
    rnd_block = (uint64_t) floor( draw_rnd(pars->rnd_gen, 0, n_blocks) );

    // And copy its content
    for(uint64_t s = 0; s < pars->boot_block_size; s++){
      block_s = block * pars->boot_block_size + s;
      rnd_block_s = rnd_block * pars->boot_block_size + s;

      if(pars->verbose >= 5)
    fprintf(stderr, "block: %lu\torig_site: %lu\trand_block:%lu\trand_site: %lu\n", block, block_s, rnd_block, rnd_block_s);

      for(uint64_t i = 0; i < pars->n_ind; i++)
    pars->geno_lkl[i][block_s] = pars->in_geno_lkl[i][rnd_block_s];
    }
  }
}

template <typename T>
struct Matrix{
  size_t x;
  size_t y;
  T** mat;
};



size_t fsize(const char* fname){
  struct stat st ;
  stat(fname,&st);
  return st.st_size;
}


Matrix<double> alloc(size_t x,size_t y){
  //fprintf(stderr,"def=%f\n",def);
  Matrix<double> ret;
  ret.x=x;
  ret.y=y;
  ret.mat= new double*[x];
  for(size_t i=0;i<x;i++)
    ret.mat[i]=new double[y];
  return ret;
}
void dalloc(Matrix<double> &ret,size_t x){
  for(size_t i=0;i<x;i++)
    delete [] ret.mat[i];
  delete [] ret.mat;
}



typedef struct emPars_t{
  int threadId; //size_t is the largest primitive datatype.
  double *inparameters;
  double *outparameters;
  Matrix<double> *GL1;
  Matrix<double> *GL2;
  int from;
  int to;
  int dim;
  double lik;
  double *sfs;//shared for all threads
  double *post;//allocated for every thread
  double *grad;//used for bfgs, length is dim-1;
} emPars;

emPars *emp = NULL;


void normalize(double *tmp,int len){
  double s=0;
  for(int i=0;i<len;i++)
    s += tmp[i];
  for(int i=0;i<len;i++)
    tmp[i] /=s;
}

double lik2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,size_t start,size_t stop){
  double res =0;
  for(size_t s=start; s<stop; s++){
    //    fprintf(stderr,"s=%d\n",s);
    double tmp =0;
    int inc =0;
    for(size_t x=0;x<GL1->y;x++)
      for(size_t y=0;y<GL2->y;y++)
    tmp += sfs[inc++]* GL1->mat[s][x]*GL2->mat[s][y];
    res +=log(tmp);
  }
  return res;
}

void emStep2(double *pre,Matrix<double> *GL1,Matrix<double> *GL2,double *post,int start,int stop,int dim){
  double inner[dim];
  for(int x=0;x<dim;x++)
    post[x] =0.0;
    
  for(int s=start; s<stop; s++){
    int inc=0;
    for(size_t x=0;x<GL1->y;x++)
      for(size_t y=0;y<GL2->y;y++){
    inner[inc] = pre[inc]*GL1->mat[s][x]*GL2->mat[s][y];
    inc++;
      }
   normalize(inner,dim);
   for(int x=0;x<dim;x++)
     post[x] += inner[x];
  }
  normalize(post,dim);
 
}


void em2(double *sfs,Matrix<double> *GL1,Matrix<double> *GL2,double tole=0.01,int maxIter=10,int dim=0){
  double oldLik,lik;
  oldLik = lik2(sfs,GL1,GL2,0,GL1->x);
  //printf("startlik=%f\n",oldLik);

  double *tmp = new double[dim];//<- wont be cleaned up, but is only allocated once
  for(int it=0; it<maxIter; it++) {
    emStep2(sfs,GL1,GL2,tmp,0,GL1->x,dim);

    for(int i=0;i<dim;i++)
      sfs[i]= tmp[i];
    
    lik = lik2(sfs,GL1,GL2,0,GL1->x);
    
    //printf("[%d] lik=%f diff=%g\n",it,lik,fabs(lik-oldLik));

    if(fabs(lik-oldLik)<tole){
      oldLik=lik;
      break;
    }
    oldLik=lik;
  }
  delete [] tmp;
}


void print(Matrix<double> &mat,FILE *fp){
  for(size_t x=0;x<mat.x;x++){
    //    fprintf(stderr,"x=%d\n",x);
    for(size_t y=0;y<mat.y;y++)
      fprintf(fp,"%f ",mat.mat[x][y]);
    fprintf(fp,"\n");
  }
}
