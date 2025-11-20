functions
{
#include spline.stan
}
data
{
  // number of stars
  int n_stars;
  // number of spline nodes
  int n_nodes;
  // number of spline nodes for the stream velocity dispersion
  int n_sig_nodes;
  // phi_1 of all the stream stars
  vector[n_stars] fi1;
  // radial velocities
  vector[n_stars] rv;
  // radial velocity errors
  vector[n_stars] erv;
  // nodes of the spline
  vector[n_nodes] fi1_nodes;
  // nodes of the spline for the velocity dispersion
  vector[n_sig_nodes] fi1_sig_nodes;
}
transformed data
{
  // the location of points in the spline nodes
  int node_ids[n_stars] = spline_findpos(fi1_nodes, fi1);
  int sig_node_ids[n_stars] = spline_findpos(fi1_sig_nodes, fi1);
}

parameters
{
  // velocities at spline nodes
  vector[n_nodes] vels;
  // log of velocity dispersions
  vector[n_sig_nodes] lsigs;
  
  // contamination fractions
  vector[n_nodes] lfracs; // these need to be transformed ex/(ex+1)
  // jacobian is also needed
  // background velocity mean
  real v0;
  // gradient in backgroud vel
  real v0grad;
  // velocity dispersion of the contamination
  real<lower=0> sig0;
}

transformed parameters
{
}

model
{

  vector[n_stars] velmod;
  vector[n_stars] lsigmod;
  vector[n_stars] lfracmod;
  vector[n_stars] lfracmod1;
  vector[n_stars] lfracmod2;
  vector[n_stars] sigmod;
  vector[n_nodes] coeffs_vel= spline_getcoeffs(fi1_nodes, vels);
  vector[n_nodes] coeffs_lfrac= spline_getcoeffs( fi1_nodes, lfracs);
  vector[n_sig_nodes] coeffs_lsig = spline_getcoeffs(fi1_sig_nodes, lsigs);
  
  // Priors
  sig0 ~ normal(110, 20);
  v0 ~ normal (0, 50);
  v0grad ~ normal(0, 0.5);
  lsigs ~ normal(1.6,1.6);
  vels ~ normal(0,300);
  
  target += sum(lfracs - 2 * log1p_exp(lfracs)); // jacobian
  
  velmod  = spline_eval(fi1_nodes,
			vels, coeffs_vel,
			fi1, node_ids);
  lfracmod  = spline_eval( fi1_nodes,
			   lfracs, coeffs_lfrac,
			   fi1, node_ids);
  lsigmod  = spline_eval(fi1_sig_nodes,
			 lsigs, coeffs_lsig,
			 fi1, sig_node_ids);
  
  sigmod = exp(lsigmod);
  lfracmod1 = lfracmod - log1p_exp(lfracmod);
  lfracmod2 = - log1p_exp(lfracmod);
  
  for (i in 1:n_stars)
    {
      target += log_sum_exp( lfracmod1[i] +
			     normal_lpdf(rv[i]|velmod[i],
					 sqrt(sigmod[i]^2 + erv[i]^2)),
			     lfracmod2[i] +
			     normal_lpdf(rv[i]|v0 + v0grad * fi1[i] ,sqrt(sig0^2+erv[i]^2)));
    }
}
