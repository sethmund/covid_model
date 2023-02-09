#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "initparams.h"
#include "helpers.h"
#include "rt_funcs.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  

// Creae funcion ha akes int .csv file and creaes
// intpus
// Name hese ex files: simulaion
// Each simulaion will have a ceraint # of runs 
// Each simulaion will have .csv files associaed 
// vv is used o deerminte he vaccintaion level
// run number is used o name he files
// Need to add dynamic vax changing



double poisson_draw(gsl_rng *r,double mu, double max_value){
    if(mu == 0){
	    return 0;
    }

    double draw_value = max_value+1;
    if(mu >= max_value){
        return max_value;
    }
    while(draw_value > max_value ){
        draw_value = gsl_ran_poisson(r,mu);
    }
    return draw_value;
}





void stoch_model(double vv, int run_number,char* fileName,struct ParameterSet p,int setting, int vax_percent){
    FILE *fptr = fopen(fileName,"w");
    double age_based_coverage[p.AGES];
    srand(time(NULL));
    gsl_rng *r;
    const gsl_rng_type *T;
    long value = rand()%10000;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r,value);

    for(int i = 0; i < p.AGES; i++){
	    p.age_based_coverage[i] = 0; 
    }
    //assigning vaccine percentages based on age
    for(int i = 0; i < p.AGES; i++){
        if ( i <= 5 & i > 0 ){
            p.age_based_coverage[i] = 0;
        }
        if ( i <= 12 & i > 5 ){
            p.age_based_coverage[i] = (.646*vv)/365.0;
        }
        if ( i <= 17 & i > 12 ){
            p.age_based_coverage[i] = (.553*vv)/365.0;
        }
        if ( i <= 49 & i > 17 ){
            p.age_based_coverage[i] = (.384*vv)/365.0;
        }
        if ( i <= 64 & i > 49 ){
            p.age_based_coverage[i] = (.506*vv)/365.0;
        }
        if ( i <= 85  & i > 64 ){
            p.age_based_coverage[i] = (.698*vv)/365.0;
        }
    }

    int vax_duration = 0;
    // ALL COMPARTMENTS
    double* m = (double*) malloc(p.AGES * sizeof(double));
    double** cm_school = (double**) malloc(p.AGES * sizeof(double*));
    double** cm_overall = (double**) malloc(p.AGES * sizeof(double*));
    double** M = (double**) malloc(p.AGES * sizeof(double*));
    double *mu =(double*) malloc(p.AGES * sizeof(double));
    double *mu_i2=(double*) malloc(p.AGES * sizeof(double));
    double* VC = (double*) malloc(p.AGES * sizeof(double));
    double* theta = (double*) malloc(p.AGES * sizeof(double));
    double* ICU_raio =(double*) malloc(p.AGES * sizeof(double)); 
    int* psi = (int*) malloc(p.years*365*sizeof(int));
    int* school_time = (int*) malloc(p.years*sizeof(int));
    double* im_prop = (double*) malloc(p.AGES * sizeof(double));
    double* S = (double*) malloc(p.AGES * sizeof(double));
    double* I1 = (double*) malloc(p.AGES * sizeof(double));
    double* VI1 = (double*) malloc(p.AGES * sizeof(double));
    double* R1 = (double*) malloc(p.AGES * sizeof(double));
    double* VR1= (double*) malloc(p.AGES * sizeof(double));
    double* VR2= (double*) malloc(p.AGES * sizeof(double));
    double* H1= (double*) malloc(p.AGES * sizeof(double));
    double* V= (double*) malloc(p.AGES * sizeof(double));
    double* I2= (double*) malloc(p.AGES * sizeof(double));
    double* VI2= (double*) malloc(p.AGES * sizeof(double*));
    double* R2= (double*) malloc(p.AGES * sizeof(double*));
    double* H2= (double*) malloc(p.AGES * sizeof(double*));
    double* Xsi1= (double*) malloc(p.AGES * sizeof(double));
    double* Xsi2= (double*) malloc(p.AGES * sizeof(double));
    double* XIV1 = (double*) malloc(p.AGES * sizeof(double));
    double* XIV2 = (double*) malloc(p.AGES * sizeof(double));
    double* XI1 = (double*) malloc(p.AGES * sizeof(double));
    double* XI2 = (double*) malloc(p.AGES * sizeof(double));
    double* XD = (double*) malloc(p.AGES * sizeof(double));
    double* N =(double*) malloc(p.AGES * sizeof(double));



    // FILENAMES
    const char *ifr_file =  "../params/ifr.csv";
    const char *vax_file =  "../params/dailyvax.csv";
    const char *m_file =  "../params/mort.csv";
    const char *n_file = "../params/us_pop.csv";
    const char *im_file = "../params/immigration_prop.csv";
    const char *overall_file = "../params/overall_contacts.csv";
    const char *icu_file = "../params/icu_ratio.csv";
    const char *school_file = "../params/school_contacts.csv";
    int psi_counter = 0;
    int perm_unvax_period = p.school_spring + p.school_break;
    int perm_vax_period = perm_unvax_period + p.perm_vax_seas_dur;
    // vaccine seasonaity loop
    for(int i = 0; i < p.years; i++){
        for(int j = 0; j < 365; j++){
            if( i == 0){
                if(j < p.vax_start){
                    psi[psi_counter] = 0; 
                }
                if(j >= p.vax_start && j < (p.vax_start + p.first_vax_seas_dur)){
                    psi[psi_counter] = 1; 
                }
                if(j >= (p.vax_start + p.first_vax_seas_dur)){
                    psi[psi_counter] = 0;
                }
            }
            if( i > 0){
                if(j < perm_unvax_period){
                    psi[psi_counter] = 0; 
                }
                if(j >= perm_unvax_period && j < perm_vax_period){
                    psi[psi_counter] = 1; 
                }
                if(j >= perm_vax_period){
                    psi[psi_counter] = 0;
                }
            }
            psi_counter += 1;
        }
    }


    float R02 = 0;
    int t = 0;
    float q1 = 0;
    float q2  = 0;
    const float q2_value = 0;
    // ALL POTENTIAL COMPARTMENTS

    double Xr1[p.AGES];
    double Y1[p.AGES];
    double Y2[p.AGES];
    double YV1[p.AGES];
    double YV2[p.AGES];
    double omega1[p.AGES];
    double omega2[p.AGES];
    double omega3[p.AGES];
    double omega4[p.AGES];
    double theta1[p.AGES];
    double theta2[p.AGES];
    double theta3[p.AGES];
    double theta4[p.AGES];
    double sv[p.AGES];
    double Xsh1[p.AGES];
    double Xsh2[p.AGES];
    double Xr1i2[p.AGES];
    double Xr2i1[p.AGES];
    double Xvvi1[p.AGES];
    double Xvvi2[p.AGES];
    double Yh1r1[p.AGES];
    double Yh2r2[p.AGES];
    double Xvr1vi2[p.AGES];
    double Xvr2vi1[p.AGES];
    double Di1[p.AGES];
    double DVi1[p.AGES];
    double DVi2[p.AGES];
    double Di2[p.AGES];
    double Dh1[p.AGES];
    double Dh2[p.AGES];
    double Vh1[p.AGES];
    double Vh2[p.AGES];
    double Vr1[p.AGES];
    double Vs[p.AGES];
    double vh1vr1[p.AGES];
    double vh2vr2[p.AGES];
    double p_vs[p.AGES];
    double r1v[p.AGES];
    double r2v[p.AGES];
    double IS[p.AGES];
    double II1[p.AGES];
    double IVI1[p.AGES];
    double IR1[p.AGES];
    double IH1[p.AGES];
    double IV[p.AGES];
    double II2[p.AGES];
    double IVI2[p.AGES];
    double IR2[p.AGES];
    double SD[p.AGES];
    double SB[p.AGES];
    double I1D[p.AGES];
    double I2D[p.AGES];
    double VI1D[p.AGES];
    double VI2D[p.AGES];
    double R2D[p.AGES];
    double R1D[p.AGES];
    double VR2D[p.AGES];
    double VR1D[p.AGES];
    double VD[p.AGES];
    double H2D[p.AGES];
    double H1D[p.AGES];
    double IH2[p.AGES];
    double D[p.AGES];
    double IVR1[p.AGES];
    double IVR2[p.AGES];
    double temp_transition = 0;
    double muIFR = 0;
    double ifr_i2_scale = 1;
    int new_yearly_imports = 100;

    // Initilizaing large datasets
    mu= initialize_unique_csv(p.AGES,ifr_file,mu);
    m = initialize_unique_csv(p.AGES,m_file,m);
    N = initialize_unique_csv(p.AGES,n_file,N);
    initialize_repeated_csv(p.AGES,vax_file,VC);
    initialize_repeated_csv(p.AGES,im_file,im_prop);
    initialize_repeated_csv(p.AGES,icu_file,ICU_raio);
    read_contact_matrices(p.AGES, overall_file,cm_overall);
    read_contact_matrices(p.AGES, school_file,cm_school);
    int counter = 0;
    int N0 = 0;
    // Setting values for theta, mu, and m
    for(int i = 0; i < p.AGES; i++){
        theta[i] = 0;
        mu[i] = (mu[i] / 1000.0)*p.IFR_mod;
        mu_i2[i] = mu[i] * ifr_i2_scale;
        m[i] = m[i];
//        m[i] = m[i] * 5.0;
        N0 += N[i];
    }
    p.N0 = N0;
    // Age-based loop for setting all transition values to zero
    for(int c = 0; c < p.AGES; c++){
        M[c] = (double*) malloc(p.AGES*sizeof(double));
        S[c] = N[c];
        for(int j = 0; j < p.AGES; j++){
            M[c][j] = cm_overall[c][j];
        }
        I1[c] = 0;
        VI1[c] = 0; 
        R1[c] = 0;
        VR1[c] = 0;
        VR2[c] = 0;
        H1[c] = 0;
        V[c] = 0;
        I2[c] = 0; 
        VI2[c] = 0;
        R2[c] = 0; 
        H2[c] = 0; 
        Xsi1[c] = 0;
        Xsi2[c] = 0;
        Xr1[c] = 0;
        Y1[c] = 0;
        Y2[c] = 0;
        YV1[c] = 0;
        YV2[c] = 0;
        omega1[c] = 0;
        omega2[c] = 0;
        omega3[c] = 0;
        omega4[c] = 0;
        theta1[c] = 0;
        theta2[c] = 0;
        theta3[c] = 0;
        theta4[c] = 0;
        sv[c] = 0;
        Xsh1[c] = 0;
        Xsh2[c] = 0;
        Xr1i2[c] = 0;
        Xr2i1[c] = 0;
        Xvvi1[c] = 0;
        Xvvi2[c] = 0;
        Yh1r1[c] = 0;
        Yh2r2[c] = 0;
        Xvr1vi2[c] = 0;
        Xvr2vi1[c] = 0;
        Di1[c] = 0;
        DVi1[c] = 0;
        DVi2[c] = 0;
        Di2[c] = 0;
        Dh1[c] = 0;
        Dh2[c] = 0;
        Vh1[c] = 0;
        Vh2[c] = 0;
        Vr1[c] = 0;
        Vs[c] = 0;
        vh1vr1[c] = 0;
        vh2vr2[c] = 0;
        p_vs[c] = 0;
        r1v[c] = 0;
        r2v[c] = 0;
        IS[c] = 0;
        II1[c] = 0;
        IVI1[c] = 0;
        IR1[c] = 0;
        IH1[c] = 0;
        IV[c] = 0;
        II2[c] = 0; 
        IVI2[c] = 0;
        IR2[c] = 0;
        SD[c] = 0;
        SB[c] = 0;
        I1D[c] = 0;
        I2D[c] = 0;
        VI1D[c] = 0;
        VI2D[c] = 0;
        R2D[c] = 0;
        R1D[c] = 0;
        VR2D[c] = 0;
        VR1D[c] = 0;
        VD[c] = 0;
        H2D[c] = 0;
        H1D[c] = 0;
        IH2[c] = 0;
        D[c] = 0;
        IVR1[c] = 0;
        IVR2[c] = 0;
    }
    t = 0;
    double total_lambda = 0;
    int rand_number = 0;
    float largest_Xsi1 = 0;
    //Time loop starts here
    while(t < p.ft){
	// School contact matrix control flow
	    if(t % (p.school_spring) == 0 ){
	        for(int i = 0; i < p.AGES; i++){
	           for(int j = 0; j < p.AGES; j++){
	    	       M[i][j] -= cm_school[i][j];
	         	}
	        }
	    }
	    if(t % (p.school_spring + p.school_break) == 0 ){
	        for(int i = 0; i < p.AGES; i++){
	           for(int j = 0; j < p.AGES; j++){
                    M[i][j] += cm_school[i][j];
	         	}
	        }
        }
	    // Variant control flow 
       // if(t == p.variant_start){
       //     R02 = p.R02;
       //      for(int i = 0 ; i < new_yearly_imports; i++){
       //          rand_number = rand() % p.AGES;
       //          if(S[rand_number] > 0){
       //              S[rand_number] = S[rand_number] -  1; 
       //              I2[rand_number] = I2[rand_number] + 1; 
       //          }
       //      }
       //      q2 = q2_value;
       //  }
           
                    for(int i = 0 ; i < new_yearly_imports; i++){
                        rand_number = rand() % p.AGES;
                        if(S[rand_number] > 0){
                            S[rand_number] = S[rand_number] - 1; 
                            I1[rand_number] = I1[rand_number] + 1; 
                        }
                    }
                   q2 = 0;  
            // Ageing Loop
        if(((t-p.school_spring) % 365 == 0)  & t != 150 ){
            S = ageing(S, p);
            I1 = ageing(I1, p);
            R1= ageing(R1, p);
            if(t > p.vax_start){
                V = ageing(V, p);
                VR1= ageing(VR1, p);
                VI1 = ageing(VI1, p);
            }
            if(t > p.variant_start){
                I2 = ageing(I2, p);
                H2= ageing(H2, p);
                VR2= ageing(VR2, p);
                VI2 = ageing(VI2, p);
            }
            // populating 0-th element of susceptibles here
            S[0] = poisson_draw(r,total(N)*p.b,total(N)); 
            int gens = 0;
            vax_duration = p.perm_vax_seas_dur;
            // Importation logic
            if(t < p.variant_start){
              for(int i = 0 ; i < new_yearly_imports; i++){
                while(gens == 0){
                  rand_number = rand() % p.AGES;
                  if(S[rand_number] > 0){
                      S[rand_number] -= 1; 
                      I1[rand_number] += 1; 
                      gens = 1;
                  }
                }
              }
            }
            else{
                 int generate = 0;
                 for(int i = 0 ; i < new_yearly_imports; i++){
                     while(generate == 0){
                         rand_number = rand() % p.AGES;

                         if(S[rand_number] > 0){
                             S[rand_number] -= 1; 
                             I1[rand_number] += 1; 
                             generate = 1;
                         }
                    }
                    generate = 0;
                 }
            }
            // Cumulatively create N
            for(int k=0; k < p.AGES; k++){
               N[k] = S[k] + I1[k] + I2[k]+ VI1[k]+ VI2[k]+ V[k]+ R1[k]+ R2[k]+ H1[k]+ H2[k]+ VR1[k]+ VR2[k];
            }
        }
        else{
            if(t  == 0){
            	for(int i = 0 ; i < new_yearly_imports; i++){
            	    rand_number = rand() % p.AGES;
            	    if(S[rand_number] > 0){
            	        S[rand_number] -= 1; 
            	        I1[rand_number] += 1; 
            	    }
            	}
                q1 = q_calc(S,I1,R1,V,N,M,mu,m,p.R01,p);
                q2 = 0;
                vax_duration = p.first_vax_seas_dur;
            }
        }

        // Stochasic Age-Transition Loop
        for(int i=0; i < p.AGES; i++){
	      // Determining attack rate!
            double lambda1 = find_lambda(q1,i,p.sigma_q1,I1,VI1,M,N,S);
           // fprintf(stderr,"Q1 is : %lf \n",q1);
           // fflush(stderr);
            total_lambda += lambda1;
            double lambda2 = find_lambda(q2,i,p.sigma_q2,I2,VI2,M,N,S);
            double lambda_mean = (S[i] * lambda1);
            Xsi1[i] =  poisson_draw(r,lambda_mean,S[i]);
            if(S[i] - Xsi1[i] >= 1){
                temp_transition = (S[i] - Xsi1[i]);
                Xsi2[i] = poisson_draw(r,temp_transition*lambda2,temp_transition);
            }
            else{
                Xsi2[i] = 0;
            }
            // S -> V compartment exit logic
            if(S[i] - Xsi1[i] - Xsi2[i] >= 1){
                temp_transition = (S[i] - Xsi1[i] - Xsi2[i]);
                sv[i] = poisson_draw(r,temp_transition*psi[t]*p.age_based_coverage[i],temp_transition);
            }
            else{
                sv[i] = 0;
            }
            // SD natural mortality exit

            if(S[i] - Xsi1[i] - Xsi2[i] - sv[i] >= 1){
                temp_transition = (S[i] - Xsi1[i] - Xsi2[i] - sv[i]);
                SD[i] = poisson_draw(r,temp_transition * m[i],temp_transition);
            }
            else{
                SD[i] = 0;
            }
           // V compartment exit logic
            Vs[i] = poisson_draw(r,V[i]*(1/p.time_of_immunity),V[i]);
            
            if(V[i] - Vs[i] >= 1){
                temp_transition = V[i] - Vs[i];
                Xvvi1[i] = poisson_draw(r,temp_transition * lambda1 * (1-p.sigma_i1),temp_transition);
            }
            else{
                Xvvi1[i] = 0;
            }
            if(V[i] - Vs[i] - Xvvi1[i] >= 1){
                temp_transition = V[i] - Vs[i] - Xvvi1[i];
                Xvvi2[i] = poisson_draw(r,temp_transition * lambda2 * (1-p.sigma_i2),temp_transition);
            }
            else{
                Xvvi2[i] = 0;
            }
            if(V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i] - DVi2[i] >= 1){
                temp_transition = V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i] - DVi2[i];
                VD[i] = poisson_draw(r,temp_transition * m[i],temp_transition);
            }
            else{
                VD[i] = 0;
            }
	        // Infected Strain One Logic

            if(I1[i] >= 1){
                Y1[i] = poisson_draw(r,I1[i]*p.gamma,I1[i]);
            }
            else{
                Y1[i] = 0;
            }

            if(I1[i] - Y1[i] >= 1){
                temp_transition = (I1[i] - Y1[i]);
                theta1[i] = poisson_draw(r,temp_transition*theta[i],temp_transition);
            }
            else{
                theta1[i] = 0;
            }

            if(I1[i] - Y1[i] - theta1[i] >= 1){
                temp_transition = (I1[i] - Y1[i] - theta1[i]);
                Di1[i] = poisson_draw(r,temp_transition*mu[i]*p.gamma,temp_transition);
            }
            else{
                Di1[i] = 0;
            }

            if(I1[i] - Y1[i] - theta1[i] - Di1[i] >= 1){
                temp_transition = (I1[i] - Y1[i] - theta1[i] - Di1[i]);
                I1D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                I1D[i] = 0;
            }

            // I2 exit logic, recovered
            if(I2[i] >= 1){
                Y2[i] = poisson_draw(r,I2[i]*p.gamma,I2[i]);
            }
            else{
                Y2[i] = 0;
            }

            // I2 exit logic, hospitalization
            if(I2[i] - Y2[i] >= 1){
                temp_transition = (I2[i] - Y2[i] );
                theta2[i] = poisson_draw(r,temp_transition*theta[i],temp_transition);
            }
            else{
                theta2[i] = 0;
            }

            if(I2[i] - Y2[i] - theta2[i] >= 1){
                temp_transition = (I2[i] - Y2[i] - theta2[i] );
                Di2[i] = poisson_draw(r,temp_transition*mu_i2[i]*p.gamma,temp_transition);
            }
            else{
                Di2[i] = 0;
            }

            if(I2[i] - Y2[i] - theta2[i] - Di2[i] >= 1 ){
                temp_transition = (I2[i] - Y2[i] - theta2[i] );
                I2D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                I2D[i] = 0;
            }

            // VI1 exit logic

            YV1[i] = poisson_draw(r,VI1[i]*p.gamma,VI1[i]);

            if(VI1[i] - YV1[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] );
                theta3[i] = poisson_draw(r,temp_transition*theta[i]*(1-p.sigma_h1),temp_transition);
            }
            else{
                theta3[i] = 0;
            }

            if(VI1[i] - YV1[i] - theta3[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] - theta3[i] );
                DVi1[i] = poisson_draw(r,temp_transition*(1-p.sigma_d1)*p.gamma*mu[i],temp_transition);
            }
            else{
                DVi1[i] = 0;
            }

            if(VI1[i] - YV1[i] - theta3[i] - DVi1[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] - theta3[i] - DVi1[i]);
                VI1D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                VI1D[i] = 0;
            }

            YV2[i] = poisson_draw(r,VI2[i]*p.gamma,VI2[i]);
            if(VI2[i] - YV2[i] >= 1){
                temp_transition = (VI2[i] - YV2[i]);
                theta4[i] = poisson_draw(r,temp_transition*theta[i]*p.sigma_h2,temp_transition);
            }
            else{
                theta4[i] = 0;
            }
            if(VI2[i] - YV2[i] - theta4[i] >= 1){
                temp_transition = (VI2[i] - YV2[i] - theta4[i]);
                DVi2[i] = poisson_draw(r,temp_transition*mu_i2[i]*p.gamma*(1-p.sigma_d2),temp_transition);
            }
            else{
                DVi2[i] = 0;
            }
            if(VI2[i] - YV2[i] - theta4[i] - DVi2[i] >= 1){
                temp_transition = (VI2[i] - YV2[i] - theta4[i]);
                VI2D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                VI2D[i] = 0;
            }

            // Recovered Exit Logic to Susceptible
            omega1[i] = poisson_draw(r,R1[i]*(1/p.time_of_waning_natural),R1[i]);

	    //Recovered One transitions
            if(R1[i] - omega1[i]>= 1){
                temp_transition = R1[i] - omega1[i];
                r1v[i] = poisson_draw(r,temp_transition*psi[t]*p.age_based_coverage[i],temp_transition);
            }
            else{
                r1v[i] = 0;
            }

            if(R1[i] - omega1[i] - r1v[i]>= 1){
                temp_transition = R1[i] - omega1[i] -r1v[i];
                Xr1i2[i] = poisson_draw(r,temp_transition*(1-p.C2)*lambda2,temp_transition);
            }
            else{
                Xr1i2[i] = 0;
            }

            if(R1[i] - omega1[i] - r1v[i] - Xr1i2[i] >= 1){
                temp_transition = R1[i] - omega1[i] -r1v[i] - Xr1i2[i];
                R1D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                R1D[i] = 0;
            }

            //R2 exit logic
            omega2[i] = poisson_draw(r,R2[i]*(1/p.time_of_waning_natural),R2[i]);
            if(R2[i] - omega2[i]>= 1){
                temp_transition = R2[i] - omega2[i];
                r2v[i] = poisson_draw(r,temp_transition*psi[t]*p.age_based_coverage[i],temp_transition);
            }
            else{
                r2v[i] = 0;
            }

            if(R2[i] - omega2[i] - r2v[i]>= 1){
                temp_transition = R2[i] - omega2[i] -r2v[i];
                Xr2i1[i] = poisson_draw(r,temp_transition*(1-p.C2)*lambda1,temp_transition);
            }
            else{
                Xr1i2[i] = 0;
            }

            if(R2[i] - omega2[i] - r2v[i] - Xr2i1[i] >= 1){
                temp_transition = R2[i] - omega2[i] -r2v[i] - Xr2i1[i];
                R2D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
            }
            else{
                R2D[i] = 0;
            }

           Xvr1vi2[i] = poisson_draw(r,VR1[i]*(1-p.VC1)*lambda2*p.sigma_i1,VR1[i]); 
     	   //Vaccine Recovered One Comaprtment Logic
           if(VR1[i] - Xvr1vi2[i] >= 1){
               temp_transition = VR1[i] - Xvr1vi2[i];
               omega3[i] = poisson_draw(r,temp_transition*(1/p.time_of_immunity),temp_transition);
           }
           else{
               omega3[i] = 0;
           }

           if(VR1[i] - Xvr1vi2[i] - omega3[i] >= 1){
               temp_transition = VR1[i] - Xvr1vi2[i] - omega3[i];
               VR1D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
           }
           else{
               VR1D[i] = 0;
           }

           Xvr2vi1[i] = poisson_draw(r,VR2[i]*(1-p.sigma_VC2)*lambda1*p.sigma_q2,VR2[i]); 
           //Vaccine Recovered Two Compartment Logic
           if(VR2[i] - Xvr2vi1[i] >= 1){
               temp_transition = VR2[i] - Xvr2vi1[i];
               omega4[i] = poisson_draw(r,temp_transition*(1/p.time_of_immunity),temp_transition);
           }
           else{
               omega4[i] = 0;
           }

           if(VR2[i] - Xvr2vi1[i] - omega4[i] >= 1){
               temp_transition = VR2[i] - Xvr2vi1[i] - omega4[i];
               VR2D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
           }
           else{
               VR2D[i] = 0;
           }
           //Hospitalization One Logic
           if(H1[i] >= 1){
               Yh1r1[i] = poisson_draw(r,H1[i]*p.zeta,H1[i]);
           }
           else{
               Yh1r1[i] = 0;
           }

           if(H1[i] - Yh1r1[i] >= 1){
               temp_transition = H1[i]-Yh1r1[i];
               H1D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
           }
           else{
               H1D[i] = 0;
           }
           //Hospitalization Two Logic

           if(H2[i] >= 1){
               temp_transition = H2[i]-Yh2r2[i];
               Yh2r2[i] = poisson_draw(r,temp_transition*p.zeta,temp_transition);
           }
           else{
               Yh2r2[i] = 0;
           }

           if(H2[i] - Yh2r2[i] >= 1){
               temp_transition = H2[i]-Yh2r2[i];
               H2D[i] = poisson_draw(r,temp_transition*m[i],temp_transition);
           }
           else{
               H2D[i] = 0;
           }
           // immigration logic

            IS[i] =  poisson_draw(r,(S[i] / N[i]) * im_prop[i],S[i]);
            IV[i] =  poisson_draw(r,(V[i] / N[i]) * im_prop[i],V[i]);
            II1[i] = poisson_draw(r,(I1[i] / N[i]) * im_prop[i],I1[i]);
            IH1[i] = poisson_draw(r,((H1[i] / N[i])) * im_prop[i],H1[i]);
            IR1[i] = poisson_draw(r,(R1[i] / N[i]) * im_prop[i],R1[i]);
            IVI1[i]= poisson_draw(r,(VI1[i] / N[i]) * im_prop[i],VI1[i]);
            if(t > p.variant_start){
                IH2[i] = poisson_draw(r,(H2[i] / N[i]) * im_prop[i],H2[i]);
                IR2[i] = poisson_draw(r,(R2[i] / N[i]) * im_prop[i],R2[i]);
                IVI2[i]= poisson_draw(r,(VI2[i] / N[i]) * im_prop[i],VI2[i]);
                IVR2[i]= poisson_draw(r,(VR2[i] / N[i]) * im_prop[i],VR2[i]);
            }

           // birth logic
            if(i == 0){
                float totalN = total(N);
                SB[i] = poisson_draw(r,totalN*p.b,totalN);
            }
        }

       if(t == 0){
           fprintf(fptr,"t,vv,age,value,sim_number,vartype\n");
        }


        for(int i = 0; i < p.AGES; i++){
            S[i] = S[i] + Vs[i] + omega1[i] + omega2[i] + omega3[i] + omega4[i] - Xsi1[i] - Xsi2[i] - sv[i] - SD[i] + SB[i] + IS[i];
            I1[i] = I1[i] + Xsi1[i] + Xr2i1[i] - Y1[i] - theta1[i] - Di1[i] - I1D[i] + II1[i];
            I2[i]  = I2[i] + Xsi2[i] + Xr1i2[i] - Y2[i] - theta2[i] - Di2[i] - I2D[i] + II2[i];
            R1[i]  = R1[i] + Y1[i] + Yh1r1[i] - r1v[i] - omega1[i] - Xr1i2[i] - R1D[i] + IR1[i];
            R2[i]  = R2[i] + Y2[i] + Yh2r2[i] - r2v[i] - omega2[i] - Xr2i1[i] - R2D[i] + IR2[i];
            VI1[i] = VI1[i] + Xvvi1[i] + Xvr2vi1[i] - theta3[i] - YV1[i] - DVi1[i] - VI1D[i] + IVI1[i];
            VI2[i] = VI2[i] + Xvvi2[i] + Xvr1vi2[i] - theta4[i] - YV2[i] - DVi2[i] - VI2D[i] + IVI2[i];
            VR1[i] = VR1[i] + YV1[i] - omega3[i] - Xvr1vi2[i] - VR1D[i] + IVR1[i];
            VR2[i] = VR2[i] + YV2[i] - omega4[i] - Xvr2vi1[i] - VR2D[i] + IVR2[i];
            V[i]   = V[i] + r1v[i] + r2v[i] + sv[i] - Xvvi1[i]  - Xvvi2[i] - Vs[i]  - VD[i] + IV[i];
            H1[i]  = H1[i] + theta1[i] + theta3[i] - Yh1r1[i] - H1D[i] + IH1[i];
            H2[i]  = H2[i] + theta2[i] + theta4[i] - Yh2r2[i] - H2D[i] + IH2[i];
            D[i]   = SD[i] + I1D[i] + I2D[i] + R1D[i] + R2D[i] + VD[i] + VI1D[i] + VI2D[i] + VR1D[i] + VR2D[i] + H1D[i] + H2D[i]; 
            N[i]   = S[i] + I1[i] + I2[i] + R1[i] + R2[i] + VI1[i] + VI2[i] + VR1[i] + VR2[i] + V[i] + H1[i] + H2[i];
            XIV1[i] = Xvvi1[i] + Xvr2vi1[i];
            XIV2[i] = Xvvi2[i] + Xvr1vi2[i];
            XI1[i] = Xsi1[i] + Xr2i1[i]; 
            XI2[i] = Xsi2[i] + Xr1i2[i];
            XD[i] = Di1[i] + Di2[i] + DVi1[i] + DVi2[i]; 

//          Age-based data-saving
        if(setting == 0){
          fprintf(fptr,"%d,%.2f,%d,%f,%d,N\n",t,vv,i,N[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,S\n",t,vv,i,S[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,I1\n",t,vv,i,I1[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,I2\n",t,vv,i,I2[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,R1\n",t,vv,i,R1[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,R2\n",t,vv,i,R2[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,V\n",t,vv,i,V[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,XIVI1\n",t,vv,i,XIV1[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,XIV2\n",t,vv,i,XIV2[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,Xsi1\n",t,vv,i,Xsi1[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,Xsi2\n",t,vv,i,Xsi2[i],run_number);
//          fprintf(fptr,"%d,%.2f,%d,%f,%d,XD\n",t,vv,i,XD[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,VR1\n",t,vv,i,VR1[i],run_number);
         // fprintf(fptr,"%d,%.2f,%d,%f,%d,VR2\n",t,vv,i,VR2[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,H1\n",t,vv,i,H1[i],run_number);
          //fprintf(fptr,"%d,%.2f,%d,%f,%d,H2\n",t,vv,i,H2[i],run_number);
        }
        // Storing less data for ages
        if(setting == 1){
          fprintf(fptr,"%d,%.2f,%d,%f,%d,N\n",t,vv,i,N[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,Xsi1\n",t,vv,i,Xsi1[i],run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,XD\n",t,vv,i,XD[i],run_number);
        }

      }

    if(setting == 2){
          fprintf(fptr,"%d,%.2f,%d,%f,%d,N\n",t,vv,90,total(N),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,S\n",t,vv,90,total(S),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,I1\n",t,vv,90,total(I1),run_number);
            fprintf(fptr,"%d,%.2f,%d,%f,%d,R1\n",t,vv,90,total(R1),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,V\n",t,vv,90,total(V),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,XIVI1\n",t,vv,90,total(XIV1),run_number);
           fprintf(fptr,"%d,%.2f,%d,%f,%d,Xsi1\n",t,vv,90,total(Xsi1),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,D\n",t,vv,90,total(D),run_number);
          fprintf(fptr,"%d,%.2f,%d,%f,%d,VR1\n",t,vv,90,total(VR1),run_number);
      }
       
        total_lambda = 0;
        // vv = dynamic_vv()
        //zero out all transitions
        for(int c = 0; c < p.AGES; c++){
	        XD[c] = 0; 
            Xsi1[c] = 0;
            Xsi2[c] = 0;
            Xr1[c] = 0;
            Y1[c] = 0;
            Y2[c] = 0;
            YV1[c] = 0;
            YV2[c] = 0;
            omega1[c] = 0;
            omega2[c] = 0;
            omega3[c] = 0;
            omega4[c] = 0;
            theta1[c] = 0;
            theta2[c] = 0;
            theta3[c] = 0;
            theta4[c] = 0;
            sv[c] = 0;
            Xsh1[c] = 0;
            Xsh2[c] = 0;
            Xr1i2[c] = 0;
            Xr2i1[c] = 0;
            Xvvi1[c] = 0;
            Xvvi2[c] = 0;
            Yh1r1[c] = 0;
            Yh2r2[c] = 0;
            Xvr1vi2[c] = 0;
            Xvr2vi1[c] = 0;
            Di1[c] = 0;
            DVi1[c] = 0;
            DVi2[c] = 0;
            Di2[c] = 0;
            Dh1[c] = 0;
            Dh2[c] = 0;
            Vh1[c] = 0;
            Vh2[c] = 0;
            Vr1[c] = 0;
            Vs[c] = 0;
            vh1vr1[c] = 0;
            vh2vr2[c] = 0;
            p_vs[c] = 0;
            r1v[c] = 0;
            r2v[c] = 0;
            IS[c] = 0;
            II1[c] = 0;
            IVI1[c] = 0;
            IR1[c] = 0;
            IH1[c] = 0;
            IV[c] = 0;
            II2[c] = 0; 
            IVI2[c] = 0;
            IR2[c] = 0;
            SD[c] = 0;
            SB[c] = 0;
            I1D[c] = 0;
            I2D[c] = 0;
            VI1D[c] = 0;
            VI2D[c] = 0;
            R2D[c] = 0;
            R1D[c] = 0;
            VR2D[c] = 0;
            VR1D[c] = 0;
            VD[c] = 0;
            H2D[c] = 0;
            H1D[c] = 0;
            IH2[c] = 0;
            D[c] = 0;
            IVR1[c] = 0;
            IVR2[c] = 0;
        }
        t += 1;
        float rt = rt_calc(S,I1,R1,V,N,M,mu,m,q1,p);
        // calculate rt-value
        int rt_age = -1;
//        vv = dynamic_vv(p.age_based_coverage,N,vax_duration,vax_percent);
//    
    }
    // now, we free all associated memory
    //free(p);
    free(m);
    free(mu_i2);
    free(VC);
    free(theta);
    free(ICU_raio); 
    free(psi);
    free(school_time);
    free(im_prop);
    free(S);
    free(I1);
    free(VI1);
    free(R1);
    free(VR1);
    free(VR2);
    free(H1);
    free(V);
    free(I2);
    free(VI2);
    free(R2);
    free(H2);
    free(Xsi1);
    free(Xsi2);
    free(XIV1);
    free(XIV2);
    free(XI1);
    free(XI2);
    free(XD);
    for(int i = 0;  i < p.AGES; i++){
        free(cm_overall[i]);
        free(cm_school[i]);
        free(M[i]);
    }
    free(cm_school);
    free(cm_overall);
    gsl_rng_free(r);
    fclose(fptr);
    return; 
}
