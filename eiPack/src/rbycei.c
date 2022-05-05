

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <Rdefines.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "eiutil.h"



SEXP
rbycei_fcn1 (
    SEXP alphamatrix,
    SEXP betaarray,
    SEXP TT,
    SEXP XX,
    SEXP tuneA,
    SEXP tuneB,
    SEXP NG,
    SEXP NP,
    SEXP Precincts,
    SEXP Lambda1,
    SEXP Lambda2,
    SEXP Sample,
    SEXP Thin,
    SEXP Burnin,
    SEXP Verbose,
    SEXP Savebeta,
    SEXP RR,
    SEXP betanames,
    SEXP XX_noisy,
    SEXP TT_noisy,
    SEXP Noise_scale,
    SEXP privacy,
    SEXP Pois_parameter
    ){

  int nProtected = 0, nComps = 0, counter = 0;

  R_len_t ii, rr, cc, tt, qq, kk; 
  SEXP lbm,ccount, hldr; 
  SEXP a_acc, b_acc, a_draws, b_draws, ccount_draws, output_list;  
  R_len_t ng, np, prec, thin, samp, burn, iters, verbose;
  double lambda_1, lambda_2;
  double aprop, acurr, asumm1, lbm_rc, aprop_ll, acurr_ll;
  double bprop, bprop_ref, bcurr, bcurr_ref, tt_ci, tt_Ci, bprop_x, bprop_ref_x, bcurr_x, bcurr_ref_x, ulim, bcurr_ll, bprop_ll;

  double prop_pop, total, curr_pop, curr_pop_ll, prop_pop_ll, noise_scale, pois_parameter;
  double prop_weighted_beta, curr_weighted_beta;
  double noise_xx, curr_xx, prop_xx;
  double noise_last_t, prop_last_t_ll, curr_last_t_ll, prop_last_t, curr_last_t, aggregate_beta, prop_total;


  ng = INTEGER(NG)[0];
  np = INTEGER(NP)[0];
  prec = INTEGER(Precincts)[0];
  lambda_1 = REAL(Lambda1)[0];
  lambda_2 = REAL(Lambda2)[0];

  samp = INTEGER(Sample)[0];
  thin = INTEGER(Thin)[0];
  burn = INTEGER(Burnin)[0];
  iters = burn + samp*thin;
  verbose = INTEGER(Verbose)[0];
  noise_scale = REAL(Noise_scale)[0];
  pois_parameter = REAL(Pois_parameter)[0];


  /*
   */

  PROTECT(hldr = allocVector(NILSXP, 1));
  ++nProtected;

  /*
  PROTECT(lbm = allocMatrix(REALSXP, ng , np));
  ++nProtected;
  */

  PROTECT(a_acc = allocVector(REALSXP, ng * np));
  ++nProtected;
  ++nComps;

  PROTECT(b_acc = allocVector(REALSXP, ng*(np-1)*prec));
  ++nProtected;
  ++nComps;


  PROTECT(ccount = allocMatrix(REALSXP, ng, np));
  ++nProtected;
  ++nComps;

  for(qq = 0; qq < ng*np; ++qq){
    REAL(a_acc)[qq] = 0; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = 0; }

  PROTECT(ccount_draws = allocMatrix(REALSXP, samp , ng*np));
  ++nProtected;
  ++nComps;


  PROTECT(a_draws = allocMatrix(REALSXP, samp , ng*np));
  ++nProtected;
  ++nComps;


  if(INTEGER(Savebeta)[0] == 0){
  PROTECT(b_draws = allocMatrix(REALSXP, samp, ng*np*prec));
  ++nProtected;
  ++nComps;
  }

  GetRNGstate();

int change_ctr = 0; int total_ctr = 0;

for(kk = 0; kk < iters; ++kk){
  for(ii = 0; ii < prec; ++ii){

    // begin new stuff
    if(REAL(privacy)[0] == 1) {

      total = 0;
      for (cc = 0; cc < np; ++cc) {
        total += REAL(TT)[cc + np*ii];
      }

      for (rr = 0; rr < ng; ++rr){

        noise_xx = rnorm(0, 0.05);
        total_ctr++;
        if (REAL(XX)[rr + ng*ii] + noise_xx >= 0 && REAL(XX)[rr + ng*ii] + noise_xx < 1) {
          
          curr_xx = REAL(XX)[rr + ng*ii];
          prop_xx = (REAL(XX)[rr + ng*ii] + noise_xx) / (1 + noise_xx);

          curr_pop = curr_xx * total;
          prop_pop = prop_xx * total;

          if (kk % 1000 == 0 && ii == 10) {
            printf("%i: Current: %f, Prop: %f\n", kk, curr_pop, prop_pop);
          }

          prop_pop_ll = // -log(lgammafn(prop_pop + 1)) + prop_pop*log(pois_parameter / ng) 
                        - (1 / (2*noise_scale*noise_scale))*pow((prop_xx * total - REAL(XX_noisy)[rr + ng*ii] * total), 2);
          curr_pop_ll = // -log(lgammafn(curr_pop + 1)) + curr_pop*log(pois_parameter / ng) 
                        - (1 / (2*noise_scale*noise_scale))*pow((curr_pop * total - REAL(XX_noisy)[rr + ng*ii] * total), 2);

          for (cc = 0; cc < np; ++cc) {
            prop_weighted_beta = REAL(betaarray)[rr + ng*cc + ng*np*ii] * prop_xx;
            curr_weighted_beta = REAL(betaarray)[rr + ng*cc + ng*np*ii] * curr_xx;
              
            prop_pop_ll += REAL(TT)[cc + ng*cc] * log(prop_weighted_beta);
            curr_pop_ll += REAL(TT)[cc + ng*cc] * log(curr_weighted_beta);
          }

          if (acc_tog(prop_pop_ll, curr_pop_ll) == 1){
            change_ctr++;
            REAL(XX)[rr + ng*ii] = prop_pop;
          }
        }
      }

      noise_last_t = rnorm(0, 5);
      prop_last_t = REAL(TT)[(np - 1) + np*ii] + noise_last_t;
      prop_total = total + noise_last_t;
      curr_last_t = REAL(TT)[(np - 1) + np*ii];

      if (prop_last_t >= 0 && prop_last_t <= total) {

        prop_last_t_ll = -(1 / (2*noise_scale*noise_scale))*pow((prop_last_t - REAL(TT_noisy)[rr + ng*ii]), 2);
        curr_last_t_ll = -(1 / (2*noise_scale*noise_scale))*pow((curr_last_t - REAL(TT_noisy)[rr + ng*ii]), 2);

        aggregate_beta = 0;
        for (rr = 0; rr < ng; ++rr) {
          aggregate_beta += REAL(XX)[rr + ng*ii] * REAL(betaarray)[rr + ng*cc + ng*np*ii];
        }

        prop_last_t_ll += log(lgammafn(prop_total + 1)) - log(lgammafn(prop_last_t + 1)) + prop_last_t * log(aggregate_beta) + (total - prop_last_t) * log(1 - aggregate_beta);
        curr_last_t_ll += log(lgammafn(total + 1)) - log(lgammafn(curr_last_t + 1)) + curr_last_t * log(aggregate_beta) + (total - curr_last_t) * log(1 - aggregate_beta);

        if (acc_tog(prop_last_t_ll, curr_last_t_ll) == 1) {
          REAL(TT)[(np - 1) + np*ii] = prop_last_t;
        }
      }

    }
    // end new stuff

    for(rr = 0; rr < ng; ++rr){
      for(cc = 0; cc < (np - 1); ++cc){

        bcurr = REAL(betaarray)[rr + ng*cc + ng*np*ii];
        bcurr_ref = REAL(betaarray)[rr + ng*(np-1) + ng*np*ii];
        ulim = bcurr + bcurr_ref;
        bprop = rnorm(bcurr, REAL(tuneB)[rr + ng*cc + ng*(np-1)*ii]);
        bprop_ref = ulim - bprop;


        if(bprop > 0 && bprop < ulim){
          tt_ci = REAL(TT)[cc + np*ii];
          tt_Ci = REAL(TT)[(np - 1) + np*ii];
          
          bcurr_x = 0; bcurr_ref_x = 0;
          for(qq = 0; qq < ng; ++qq){
            bcurr_x += REAL(betaarray)[qq + ng*cc + ng*np*ii]* REAL(XX)[qq + ng*ii];
            bcurr_ref_x += REAL(betaarray)[qq + ng*(np-1) + ng*np*ii] * REAL(XX)[qq + ng*ii];
          }
         

          bprop_x = bcurr_x - bcurr*REAL(XX)[rr + ng*ii] + bprop*REAL(XX)[rr + ng*ii];
          bprop_ref_x =  bcurr_ref_x - bcurr_ref*REAL(XX)[rr + ng*ii] + bprop_ref*REAL(XX)[rr + ng*ii];

          bprop_ll = beta_ll(bprop, bprop_ref, REAL(alphamatrix)[rr + ng*cc], REAL(alphamatrix)[rr + ng*(np-1)],tt_ci, tt_Ci, bprop_x, bprop_ref_x);
          bcurr_ll = beta_ll(bcurr, bcurr_ref,  REAL(alphamatrix)[rr + ng*cc], REAL(alphamatrix)[rr + ng*(np-1)], tt_ci, tt_Ci, bcurr_x, bcurr_ref_x);
        
          /* Rprintf("%f %f %f %f %f\n", bprop, bprop_ref, bprop_x, bprop_ref_x, bprop_ll - bcurr_ll); 
          */

          if(acc_tog(bprop_ll, bcurr_ll) == 1){
            REAL(betaarray)[rr + ng*cc + ng*np*ii] = bprop;
            REAL(betaarray)[rr + ng*(np-1) + ng*np*ii] = bprop_ref;
            REAL(b_acc)[rr + ng*cc + ng*(np - 1)*ii] += 1;
        
          }
        }
      }
    }
  }


  PROTECT(lbm = logbetamat(betaarray, NG, NP, Precincts));

  for(rr = 0; rr < ng; ++rr){
    for(cc = 0; cc < np; ++cc){
      
      acurr = REAL(alphamatrix)[rr + ng*cc];
      aprop = rnorm(acurr, REAL(tuneA)[rr + ng*cc]);
     
      lbm_rc = REAL(lbm)[rr + ng*cc];
      asumm1 = 0; 
      for(tt = 0; tt < np; ++tt){
        asumm1 += REAL(alphamatrix)[rr + ng*tt];
      }
      asumm1 = asumm1 - acurr;

      if(aprop > 0){
        aprop_ll = alpha_ll(aprop, asumm1, lbm_rc, lambda_1, lambda_2, prec);
        acurr_ll = alpha_ll(acurr, asumm1, lbm_rc, lambda_1, lambda_2, prec);
        /*Rprintf("%f %f \n", aprop, aprop_ll - acurr_ll);
         */
        if(acc_tog(aprop_ll, acurr_ll) == 1){
          REAL(alphamatrix)[rr + ng*cc] = aprop;
          REAL(a_acc)[rr + ng*cc] += 1;
        }
      }
    }
  }
  UNPROTECT(1);



  if(kk >= burn && ((kk % thin) == 0)){


    ccount = cellcount(betaarray,RR, NG, NP, Precincts);


    for(qq = 0; qq < np*ng; ++qq){
      REAL(a_draws)[counter + qq*samp] = REAL(alphamatrix)[qq];
      REAL(ccount_draws)[counter + qq*samp] = REAL(ccount)[qq];
    }

    if(INTEGER(Savebeta)[0] == 0){
      for(qq = 0; qq < np*ng*prec; ++qq){
       REAL(b_draws)[counter + qq*samp] = REAL(betaarray)[qq];
      }
    }


    if(INTEGER(Savebeta)[0] == 2){
      write_beta(betaarray, betanames);
    }

     counter += 1;
  }

  if(verbose > 0 && kk % verbose == 0){
    Rprintf("\n MCMC iteration %i of %i \n", kk + 1, iters); 
  }

  R_CheckUserInterrupt();

}

   /*
    *   PROTECT(dim_matrix1 = allocVector(INTSXP, 2));
    *++nProtected; 
    *   INTEGER(dim_matrix1)[0] = ncol;
    *   INTEGER(dim_matrix1)[1] = nrow;
    *   
    *   setAttrib(matrix1, R_DimSymbol, dim_matrix1);
    */


for(qq = 0; qq < ng*np; ++qq){
    REAL(a_acc)[qq] = REAL(a_acc)[qq]/iters; }
  for(qq = 0; qq < ng*(np-1)*prec;++qq){
    REAL(b_acc)[qq] = REAL(b_acc)[qq]/iters; }


  if(INTEGER(Savebeta)[0]==0){
    PROTECT(output_list = allocVector(VECSXP, 5));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, a_draws);
    SET_VECTOR_ELT(output_list, 1, b_draws);
    SET_VECTOR_ELT(output_list, 2, a_acc);    
    SET_VECTOR_ELT(output_list, 3, b_acc);
    SET_VECTOR_ELT(output_list, 4, ccount_draws);
  }else{
    PROTECT(output_list = allocVector(VECSXP, 4));
    ++nProtected; 
    SET_VECTOR_ELT(output_list, 0, a_draws);
    SET_VECTOR_ELT(output_list, 1, a_acc);    
    SET_VECTOR_ELT(output_list, 2, b_acc);
        SET_VECTOR_ELT(output_list, 3, ccount_draws);
  }


  PutRNGstate();
  UNPROTECT(nProtected); 

  // printf("%i, %i\n", change_ctr, total_ctr);

  return(output_list); 

}
