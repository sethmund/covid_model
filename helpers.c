#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simconstants.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  
#include "initparams.h"
#include "helpers.h"

//sigma = reduction of infectiousness
double find_lambda(double q,int age, double sigma, double* I, double* VI, double** M, double* N,double* S){
    double row_sum = 0;
    double infec_contact = 0; 
    double breakthrough_contact = 0;
   //         fprintf(stderr,"VI at 0: %lf \n",VI[0]);
   //         fflush(stderr);
     for(int i = 0; i < AGES; i++){
         infec_contact += (M[age][i]*(I[i]/N[i]));
     }
     for(int i = 0; i < AGES; i++){
         breakthrough_contact += (M[age][i]*(VI[i]/N[i]));
     }
     row_sum = q*(infec_contact + sigma*breakthrough_contact);
    return row_sum;
}


//Creating vaccination array 
//double* create_vaccine(int* psi,int vv,){
//    double age_based_coverage[AGES];
//    for(int i = 0; i < AGES; i++){
//	    age_based_coverage[i] = 0; 
//    }
//    //assigning vaccine percentages based on age
//    for(int i = 0; i < AGES; i++){
//        fprintf(stderr,"VV VALUE: %f \n",vv);
//        fflush(stderr);
//        if ( i <= 5 & i > 0 ){
//            age_based_coverage[i] = 0;
//        }
//        if ( i <= 12 & i > 5 ){
//            age_based_coverage[i] = (.646*vv)/365.0;
//        }
//        if ( i <= 17 & i > 12 ){
//            age_based_coverage[i] = (.553*vv)/365.0;
//        }
//        if ( i <= 49 & i > 17 ){
//            age_based_coverage[i] = (.384*vv)/365.0;
//        }
//        if ( i <= 64 & i > 49 ){
//            age_based_coverage[i] = (.506*vv)/365.0;
//        }
//        if ( i <= 85  & i > 64 ){
//            age_based_coverage[i] = (.698*vv)/365.0;
//        }
//    }
//    return age_based_coverage;
//}
//int* set_psi(int vaccine_start,int first_vsd,int perm_up,int perm_vp){
//    //setting psi value
//    for(int i = 0; i < years; i++){
//        for(int j = 0; j < 365; j++){
//            if( i == 0){
//                if(j < vaccine_start){
//                    psi[psi_counter] = 0; 
//                }
//                if(j >= vaccine_start && j < (vaccine_start + first_vsd)){
//                    psi[psi_counter] = 1; 
//                }
//                if(j >= (vaccine_start + first_vsd)){
//                    psi[psi_counter] = 0;
//                }
//            }
//            if( i > 0){
//                if(j < perm_up){
//                    psi[psi_counter] = 0; 
//                }
//                if(j >= perm_up && j < perm_vp){
//                    psi[psi_counter] = 1; 
//                }
//                if(j >= perm_vp){
//                    psi[psi_counter] = 0;
//                }
//            }
//            psi_counter += 1;
//        }
//    }
//
//
//} 

double* ageing(double* L,struct ParameterSet p){
    double *new_L = (double*) malloc(AGES*sizeof(double));
    // this will get set in the model function
    new_L[0] = 0;
    for(int i = 1; i < p.AGES; i++){
        if(i == 84){
            new_L[i] = L[i-1] + L[i];
        }
        else{
            new_L[i] = L[i-1];
        }
    }
    free(L);
    return new_L;
}

double total(double* L){
    double sum = 0;
    for(int i = 0; i<AGES; i++){
        sum += L[i];
    }
    return sum;
}

//double poisson_draw(double mu, double max_value){
//
//    double draw_value = max_value+1;
//    // gsl_rng *test; 
//    // const gsl_rng_type *T;
//    // gsl_rng_env_setup();
//    // T = gsl_rng_default;
//    // test = gsl_rng_alloc(T);
//    if(mu >= max_value){
//        return max_value;
//    }
//    while(draw_value > max_value ){
//        draw_value = gsl_ran_poisson(r,mu);
//    }
//    return draw_value;
//}



/*
Dynamic_vv calculates a new vv value based on the population at the time, the vaccine regime
the desired vaccine coverage percentage 
*/
double dynamic_vv(double* ABR, double* N, int VR, int ideal_vax_coverage){
    float total_N = total(N);
    int counter = 1; 
    double vv_val= 0;
    float test_WVC = 0;
    float vax_decimal = ideal_vax_coverage/100.0;
    for(int i =100; i < 13000; i++){
        double test_vv = i / 1000;
        for(int j = 0; j < AGES; j++){
            test_WVC += ABR[j]*N[j];
        }
        test_WVC = (test_WVC*test_vv*VR);
        test_WVC = test_WVC / total(N);
            fprintf(stderr,"TEST WVC: %lf \n",test_WVC);
            fflush(stderr);
        double diff = abs(vax_decimal - test_WVC);
        fprintf(stderr,"TEST WVC: %lf, VAX DECIMAL: %lf \n",test_WVC,vax_decimal);
        fflush(stderr);
        if(diff < .0001){
            vv_val = test_vv;
            counter += 1;
            fprintf(stderr,"VV done %lf !\n",vv_val);
            fflush(stderr);
            break;
        }
    }
    return vv_val;
}
