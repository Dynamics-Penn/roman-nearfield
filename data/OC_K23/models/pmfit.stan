functions
{
#include spline.stan
}
data
{
  // number of stars
  int n_stars;
  // number of spline nodes
  int n_pm_nodes;
  // number of spline nodes
  int n_frac_nodes;
  int n_lsig_nodes;
  // phi_1 of all the stream stars
  vector[n_stars] fi1;

  // radial velocities
  vector[n_stars] pm;
  // radial velocity errors
  vector[n_stars] epm;

  // nodes of the spline
  vector[n_pm_nodes] pm_nodes;
  vector[n_frac_nodes] frac_nodes;
  vector[n_lsig_nodes] lsig_nodes;
  vector[n_stars] pmcut1;
  vector[n_stars] pmcut2;

}
transformed data
{
  // the location of points in the spline nodes
  int node_pm_ids[n_stars] = spline_findpos(pm_nodes, fi1);
  int node_frac_ids[n_stars] = spline_findpos(frac_nodes, fi1);
  int node_lsig_ids[n_stars] = spline_findpos(lsig_nodes, fi1);
}

parameters
{
  // velocities at spline nodes
  vector[n_pm_nodes] pms;

  // contamination fractions
  vector[n_frac_nodes] lfracs; // these need to be transformed 0.5 + arctan(x)/pi
  vector[n_lsig_nodes] lsigs; // log sigma
  // jacobian is also needed
  // background velocity mean
}

transformed parameters
{
}

model
{

  vector[n_stars] pmmod;
  vector[n_stars] lfracmod;
  vector[n_stars] fracmod;
  vector[n_stars] sigmod;
  vector[n_stars] fractmp;
  vector[n_stars] lfracmod1;
  vector[n_stars] lfracmod2;
  vector[n_pm_nodes] coeffs_pm = spline_getcoeffs(pm_nodes, pms);
  vector[n_frac_nodes] coeffs_lfrac = spline_getcoeffs(frac_nodes, lfracs);
  vector[n_lsig_nodes] coeffs_lsig = spline_getcoeffs(lsig_nodes, lsigs);

  // Priors
  pms ~ normal(0, 7);
  //lsigs ~ normal(-3, 3);
  lsigs ~ normal(-2.5,1.5);

  // the jacobian
  for (i in 1:n_frac_nodes)
    {
      target += lfracs[i] - 2 * log_sum_exp(lfracs[i], 0);

    }
  pmmod  = spline_eval(pm_nodes, pms, coeffs_pm,
	              fi1, node_pm_ids);
  lfracmod  = spline_eval(frac_nodes, lfracs, coeffs_lfrac,
			          fi1, node_frac_ids);
  sigmod  = exp(spline_eval(lsig_nodes, lsigs, coeffs_lsig,
              	              fi1, node_lsig_ids));
  for (i in 1:n_stars)
    {
      fractmp[i]  = log_sum_exp(lfracmod[i], 0);
    } // this is log(e^x+1)
  fracmod = exp(lfracmod - fractmp);
  lfracmod1 = lfracmod - fractmp;
  lfracmod2 =  - fractmp;
  for (i in 1:n_stars)
    {
      real cursig = sqrt(sigmod[i] ^2 + epm[i]^2);
      real lossfrac = normal_cdf(pmcut2[i]| pmmod[i], cursig) -
       normal_cdf(pmcut1[i]|pmmod[i], cursig);
      real renorm = lossfrac * fracmod[i] + (1 - fracmod[i]);
      target += log_sum_exp(lfracmod1[i] +
			          normal_lpdf(pm[i]|pmmod[i], cursig),

			    lfracmod2[i] + log(1/(pmcut2[i]-pmcut1[i]))
			    ) - log(renorm);
    }
   //print(target(), lfracs, lsigs);
}
