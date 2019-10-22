// SIS-based evidence synthesis for M. genitalium in the UK, with PID compartment

data {

    int<lower=0> prev_num_ct;      // numerator, for prevalence data
    int<lower=0> prev_denom_ct;    // denominator, for prevalence data

    int<lower=0> prev_num_mg;      // numerator, for prevalence data
    int<lower=0> prev_denom_mg;    // denominator, for prevalence data
    
// Only MG has persistence data. CT persistence time is done with an informative prior.
    int<lower=0> N_persist; // number of data points for persistence time
    real<lower=0> persist_t; // times for persistence follow-up
    int<lower=0> persist_num[N_persist]; // numerators for persistence follow-up
    int<lower=0> persist_denom[N_persist]; // denominators for persistence follow-up
    
    int<lower=0> pid_given_infected_num_ct; // number of M gen-infected women in cohort study who developed PID over one year follow-up
    int<lower=0> pid_given_infected_denom_ct; // number of M gen-infected women in cohort study
    int<lower=0> pid_given_susceptible_num_ct; // number of M gen-uninfected women in cohort study who developed PID over one year follow-up
    int<lower=0> pid_given_susceptible_denom_ct; // number of M gen-uninfected women in cohort study
  
    int<lower=0> pid_given_infected_num_mg; // number of M gen-infected women in cohort study who developed PID over one year follow-up
    int<lower=0> pid_given_infected_denom_mg; // number of M gen-infected women in cohort study
    int<lower=0> pid_given_susceptible_num_mg; // number of M gen-uninfected women in cohort study who developed PID over one year follow-up
    int<lower=0> pid_given_susceptible_denom_mg; // number of M gen-uninfected women in cohort study
  
}

parameters{

  real<lower=0> alpha_SC; // FoI Ct
  real<lower=0> alpha_SM; // FoI Mgen
  
  real<lower=0> alpha_CS; // recovery rate Ct
  real<lower=0> alpha_MS; // recovery rate Mgen
  
  real<lower=0> alpha_SP; // PID progression - "background" rate
  real<lower=0> alpha_CP; // PID progression - Ct
  real<lower=0> alpha_MP; // PID progression - Mgen

}

transformed parameters {

  real C_star; // steady-state Ct prevalence
  real M_star; // steady-state Mgen prevalence
  
  real p_persist_M; // probability of persistence - Mgen

  matrix[5,5] Q;
  
  matrix[5,5] P_1yr; // one-year transition matrix
  
  C_star = alpha_SC / (alpha_CS + alpha_SC);
  M_star = alpha_SM / (alpha_MS + alpha_SM);
  
  p_persist_M = M_star + (1 - M_star)*exp(-(alpha_MS + alpha_SM)*persist_t);

  Q = [[-(alpha_SC + alpha_SM + alpha_SP), alpha_SC, alpha_SM, 0, alpha_SP],
        [alpha_CS, -(alpha_CS + alpha_SM + alpha_SP + alpha_CP), 0, alpha_SM, alpha_SP + alpha_CP],
        [alpha_MS, 0, -(alpha_MS + alpha_SC + alpha_SP + alpha_MP), alpha_SC, alpha_SP + alpha_MP],
        [0, alpha_MS, alpha_CS, -(alpha_MS + alpha_CS + alpha_SP + alpha_CP + alpha_MP), alpha_SP + alpha_CP + alpha_MP],
        [0, 0, 0, 0, 0]];

  P_1yr = matrix_exp(Q);   

}

model {

  // priors

  alpha_SC ~ gamma(1, 2); // FoI Ct
  alpha_SM ~ gamma(1, 2); // FoI Mgen
  
  alpha_CS ~ normal(0.74, (0.89-0.61)/3.919928); // recovery rate Ct
  alpha_MS ~ gamma(1, 2); // recovery rate Mgen
  
  alpha_SP ~ gamma(1, 2); // PID progression - "background" rate
  alpha_CP ~ gamma(1, 2); // PID progression - Ct
  alpha_MP ~ gamma(1, 2); // PID progression - Mgen

  // likelihood
  persist_num ~ binomial(persist_denom, p_persist_M);

  prev_num_ct ~ binomial(prev_denom_ct, C_star);
  prev_num_mg ~ binomial(prev_denom_mg, M_star);
    
  pid_given_infected_num_ct ~ binomial(pid_given_infected_denom_ct, (1 - M_star)*P_1yr[2,5] + M_star*P_1yr[4,5]); // NB infections are independent, so can use M_star and (1 - M_star)
  pid_given_infected_num_mg ~ binomial(pid_given_infected_denom_mg, (1 - C_star)*P_1yr[3,5] + C_star*P_1yr[4,5]); // NB infections are independent, so can use C_star and (1 - C_star)
    
  pid_given_susceptible_num_ct ~ binomial(pid_given_susceptible_denom_ct, (1 - M_star)*P_1yr[1,5] + M_star*P_1yr[3,5]);
  pid_given_susceptible_num_mg ~ binomial(pid_given_susceptible_denom_mg, (1 - C_star)*P_1yr[1,5] + C_star*P_1yr[2,5]);
    
}

generated quantities {

    real<lower=0,upper=1> natsal_prev_mg[4]; // four age groups
    real<lower=0,upper=1> natsal_prev_ct[4]; 
    
    real<lower=0> natsal_pid[4];
    real<lower=0> natsal_pid_obs[4]; // account for probability that case recorded in surveillance
    real<lower=0> natsal_pid_mg[4];
    real<lower=0> natsal_pid_mg_obs[4]; // account for probability that case recorded in surveillance
    real<lower=0> natsal_pid_ct[4];
    real<lower=0> natsal_pid_ct_obs[4]; // account for probability that case recorded in surveillance

    natsal_prev_mg[1] = beta_rng(5+1, 211-5+1); // based on Sonnenberg IJE 2015
    natsal_prev_mg[2] = beta_rng(5+1, 379-5+1); 
    natsal_prev_mg[3] = beta_rng(11+1, 799-11+1); 
    natsal_prev_mg[4] = beta_rng(9+1, 859-9+1); 
    
    natsal_prev_ct[1] = beta_rng(3.667971+1, 101.313295+1); // based on Natsal-2
    natsal_prev_ct[2] = beta_rng(7.104704+1, 252.014263+1); 
    natsal_prev_ct[3] = beta_rng(11.782488+1, 677.866813+1); 
    natsal_prev_ct[4] = beta_rng(4.046203+1, 667.071823+1); 
    
    for(i in 1:4){
      // pid from all causes
      natsal_pid[i] = 100000*( 
        (1-natsal_prev_mg[i])*(1-natsal_prev_ct[i]) * alpha_SP + // susceptible
        natsal_prev_mg[i] * (1-natsal_prev_ct[i]) * (alpha_SP + alpha_MP) + // Mg only
        natsal_prev_ct[i] * (1-natsal_prev_mg[i]) * (alpha_SP + alpha_CP) + // Ct only
        natsal_prev_mg[i] * natsal_prev_ct[i] * (alpha_SP + alpha_MP + alpha_CP) // both infections
        ); 
        
      natsal_pid_obs[i] = natsal_pid[i] * beta_rng(12 + 9 + 1, 15 + 23 - 12 - 9 + 1);
        
      // pid attributable to Mg
      natsal_pid_mg[i] = 100000*natsal_prev_mg[i] * (alpha_MP);
      natsal_pid_mg_obs[i] = natsal_pid_mg[i] * beta_rng(12 + 9 + 1, 15 + 23 - 12 - 9 + 1);
      
      // pid attributable to Ct
      natsal_pid_ct[i] = 100000*natsal_prev_ct[i] * (alpha_CP);
      natsal_pid_ct_obs[i] = natsal_pid_ct[i] * beta_rng(12 + 9 + 1, 15 + 23 - 12 - 9 + 1);
      }

}