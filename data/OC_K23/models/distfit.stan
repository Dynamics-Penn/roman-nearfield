functions
{
#include spline.stan
}
data 
{
  // number of stars
  int n_stars;
  // phi_1 of all the stream stars
  vector[n_stars] fi1;
  // radial velocities
  vector[n_stars] dm;
  // radial velocity errors
  vector[n_stars] edm;
  vector[n_stars] dmcut1;
  vector[n_stars] dmcut2;

  // number of spline nodes 
  int n_nodes;
  // number of sigma nodes
  int n_sig_nodes;
  // nodes of the spline
  vector[n_nodes] fi1_nodes;
  // nodes of the spline
  vector[n_sig_nodes] sig_nodes;
  
}
transformed data
{
  // the location of points in the spline nodes
  int node_ids[n_stars] = spline_findpos(fi1_nodes, fi1);
  int node_sig_ids[n_stars] = spline_findpos(sig_nodes, fi1);
}

parameters
{
  // velocities at spline nodes
  vector[n_nodes] dms;
  // velocities at spline nodes
  vector[n_sig_nodes] lsigs;
  // contamination fractions
  vector[n_nodes] lfracs; // these need to be transformed ex/(ex+1)
  // jacobian is also needed
  // background velocity mean
}

transformed parameters
{
}

model
{

  vector[n_stars] dmmod;
  vector[n_stars] lsigmod;
  vector[n_stars] sigmod;
  vector[n_stars] lfracmod;
  vector[n_stars] lfracmod1;
  vector[n_stars] lfracmod2;

  vector[n_nodes] coeffs_dm = spline_getcoeffs(fi1_nodes, dms);
  vector[n_nodes] coeffs_lfrac = spline_getcoeffs(fi1_nodes, lfracs);
  vector[n_sig_nodes] coeffs_sig = spline_getcoeffs(sig_nodes, lsigs);
  real xsig;
  
  // Priors
  dms ~ normal(18, 2);
  lsigs ~ normal(log(0.1), 2);
  
  target += sum(lfracs-2*log1p_exp(lfracs)); // jacobian
  
  dmmod  = spline_eval(fi1_nodes,
		       dms, coeffs_dm,
		       fi1, node_ids);
  lfracmod  = spline_eval(fi1_nodes,
			  lfracs, coeffs_lfrac,
			  fi1, node_ids);
  lsigmod  = spline_eval(sig_nodes,
			 lsigs, coeffs_sig,
			 fi1, node_sig_ids);
  sigmod = exp(lsigmod);
  lfracmod1 = lfracmod - log1p_exp(lfracmod);
  lfracmod2 = - log1p_exp(lfracmod);
  
  for (i in 1:n_stars)
    {
      xsig = sqrt(sigmod[i]*sigmod[i] + edm[i]*edm[i]);
      target += log_sum_exp( lfracmod1[i] +
			     normal_lpdf(dm[i]|dmmod[i], xsig)-
			     log(1e-100+
				 (normal_cdf(dmcut2[i],dmmod[i], xsig)-
				  normal_cdf(dmcut1[i],dmmod[i], xsig))),
			     lfracmod2[i] + log(1/(dmcut2[i]-dmcut1[i]))
			     );
    }
  //    print(target());
}
