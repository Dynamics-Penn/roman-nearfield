functions
{
#include spline_precompute.stan

  real partial_lpdf_sum(int[] hh_slice,
			int start,
			int end,
			
			vector int_nodes,
                        vector ints,
			vector coeffs_int,
                        matrix mat_int,
			int[] node_ids_int,
			
			vector fi2_nodes,
                        vector fi2s,
			vector coeffs_fi2,
                        matrix mat_fi2,
			int[] node_ids_fi2,
			
			vector logwidth_nodes,
                        vector logwidths,
			vector coeffs_logwidth,
                        matrix mat_logwidth,
			int[] node_ids_logwidth,
			
			vector bg_nodes,
                        vector bgs,
			vector coeffs_bg,
                        matrix mat_bg,
			int[] node_ids_bg,
			
			vector bgsl_nodes,
                        vector bgsls,
			vector coeffs_bgsl,
                        matrix mat_bgsl,
			int[] node_ids_bgsl,
			
			vector bgsl2_nodes,
                        vector bgsl2s,
			vector coeffs_bgsl2,
                        matrix mat_bgsl2,
			int[] node_ids_bgsl2,
			
                        vector bgpm1s,
			vector coeffs_bgpm1,
                        vector bgslpm1s,
			vector coeffs_bgslpm1,
			
                        vector bgsl2pm1s,
			vector coeffs_bgsl2pm1,
			
                        vector bgpm2s,
			vector coeffs_bgpm2,
			
                        vector bgslpm2s,
			vector coeffs_bgslpm2,
			
                        vector bgsl2pm2s,
			vector coeffs_bgsl2pm2,
			
			vector pm_nodes,
                        vector pm1s,
			vector coeffs_pm1,
                        matrix mat_pm,
			int[] node_ids_pm,
			
                        vector pm2s,
			vector coeffs_pm2,
			
			vector pmlogwidth_nodes,
                        vector pmlogwidths,
			vector coeffs_pmlogwidth,
                        matrix mat_pmlogwidth,
			int[] node_ids_pmlogwidth,
			
			vector x,
			vector y,
			vector dpm1,
			vector dpm2
			)
			
  {
    int ncur = end-start+1;
    vector[ncur] int_val;
    vector[ncur] logwidth_val;
    vector[ncur] bg_val;
    vector[ncur] bgsl_val;
    vector[ncur] bgsl2_val;
    vector[ncur] bgpm1_val;
    vector[ncur] bgslpm1_val;
    vector[ncur] bgsl2pm1_val;
    vector[ncur] bgpm2_val;
    vector[ncur] bgslpm2_val;
    vector[ncur] bgsl2pm2_val;
    vector[ncur] fi2_val;
    vector[ncur] pm1_val;
    vector[ncur] pm2_val;
    vector[ncur] pmlogwidth_val;
    vector[ncur] xmod;
    vector[ncur] logbg_pix;
    vector[ncur] int_pix;
    

    // 6 locations for var names (in RHS)
    // Actual evaluation of the splines at each pixel
    int_val = spline_eval(int_nodes,
			  ints, coeffs_int,
			  mat_int[start:end], node_ids_int[start:end]);
    
    fi2_val = spline_eval(fi2_nodes,
			  fi2s, coeffs_fi2,
			  mat_fi2[start:end], node_ids_fi2[start:end]);
    
    logwidth_val = spline_eval(logwidth_nodes,
			       logwidths, coeffs_logwidth,
			       mat_logwidth[start:end], node_ids_logwidth[start:end]);
    
    bg_val = spline_eval(bg_nodes,
			 bgs, coeffs_bg,
			 mat_bg[start:end], node_ids_bg[start:end]);
    
    bgsl_val = spline_eval(bgsl_nodes,
			   bgsls, coeffs_bgsl,
			   mat_bgsl[start:end], node_ids_bgsl[start:end]);
    
    bgsl2_val = spline_eval(bgsl2_nodes,
			    bgsl2s, coeffs_bgsl2,
			    mat_bgsl2[start:end], node_ids_bgsl2[start:end]);
    
    // here we are reusing mat_ ids and nodes so
    // we need to provide different coeffs and node values
    bgpm1_val = spline_eval(bg_nodes,
			    bgpm1s, coeffs_bgpm1,
			    mat_bg[start:end], node_ids_bg[start:end]);
    
    bgslpm1_val = spline_eval(bgsl_nodes,
			      bgslpm1s, coeffs_bgslpm1,
			      mat_bgsl[start:end], node_ids_bgsl[start:end]);
    
    bgsl2pm1_val = spline_eval(bgsl2_nodes,
			       bgsl2pm1s, coeffs_bgsl2pm1,
			       mat_bgsl2[start:end], node_ids_bgsl2[start:end]);
    
    bgpm2_val = spline_eval(bg_nodes,
			    bgpm2s, coeffs_bgpm2,
			    mat_bg[start:end], node_ids_bg[start:end]);
    
    bgslpm2_val = spline_eval(bgsl_nodes,
			      bgslpm2s, coeffs_bgslpm2,
			      mat_bgsl[start:end], node_ids_bgsl[start:end]);
    
    bgsl2pm2_val = spline_eval(bgsl2_nodes,
			       bgsl2pm2s, coeffs_bgsl2pm2,
			       mat_bgsl2[start:end], node_ids_bgsl2[start:end]);
    
    pm1_val = spline_eval(pm_nodes,
			  pm1s, coeffs_pm1,
			  mat_pm[start:end], node_ids_pm[start:end]);
    
    pm2_val = spline_eval(pm_nodes,
			  pm2s, coeffs_pm2,
			  mat_pm[start:end], node_ids_pm[start:end]);
    
    pmlogwidth_val = spline_eval(pmlogwidth_nodes,
				 pmlogwidths, coeffs_pmlogwidth,
				 mat_pmlogwidth[start:end], node_ids_pmlogwidth[start:end]);

    // log densities of the background/stream at each pixel
    logbg_pix = (
		 (bg_val + bgsl_val .* (y[start:end]/10) +
		  bgsl2_val .* square(y[start:end]/10)) +
		 (bgpm1_val + bgslpm1_val .* (y[start:end]/10) +
		  bgsl2pm1_val .* square(y[start:end]/10)) .* dpm1[start:end] +
		 (bgpm2_val + bgslpm2_val .* (y[start:end]/10) +
		  bgsl2pm2_val .* square(y[start:end]/10)) .* dpm2[start:end]);
    int_pix = (int_val
	       - 0.5 * square(y[start:end]-fi2_val) ./ exp(2 * logwidth_val)
	       - 0.5 * (square(dpm1[start:end] - pm1_val) +
			square(dpm2[start:end] - pm2_val))
	       ./ exp(2 * pmlogwidth_val));

    for (i in 1:(end-start+1))
    {
      xmod[i] = log_sum_exp(logbg_pix[i], int_pix[i]);
    }

    //Likelihood
    return poisson_log_lpmf(hh_slice|xmod);
  
  }  
}

data {
  // number of pixels in the density map
  int n_pix;
  // number counts of objects in pixels
  int hh[n_pix];
  // the x locations of pixels
  vector[n_pix] x;
  // the y locations of pixels
  vector[n_pix] y;
  vector[n_pix] dpm1;
  vector[n_pix] dpm2;
  
  
  // number of nodes in the intencity spline
  int n_int_nodes;
  // location of nodes or the intensity spline
  vector[n_int_nodes] int_nodes;
  
  // number of nodes in the spatial density spline
  int n_fi2_nodes;
  // location of nodes for the spatial density spline
  vector[n_fi2_nodes] fi2_nodes;
  
  // number of nodes in the stream width spline
  int n_logwidth_nodes;
  // location of nodes for the width spline
  vector[n_logwidth_nodes] logwidth_nodes;
  
  // number of nodes for the background density spline
  int n_bg_nodes;
  // location of nodes the background density spline
  vector[n_bg_nodes] bg_nodes;

  // the same thing as above for the background slope
  int n_bgsl_nodes;
  vector[n_bgsl_nodes] bgsl_nodes;
  
  // the same thing as above for the quadratic background slope
  int n_bgsl2_nodes;
  vector[n_bgsl2_nodes] bgsl2_nodes;
  
  int n_pm_nodes;
  vector[n_pm_nodes] pm_nodes;

  int n_pmlogwidth_nodes;
  vector[n_pmlogwidth_nodes] pmlogwidth_nodes;
  int grainsize;

}

transformed data
{
  // These are the calculation of in which node interal each pixel
  // should go to. Since the location of nodes is different for the 
  // density/track/width etc the locations ids will be differenet for
  // each of those parameters

  int node_ids_int[n_pix] = spline_findpos(int_nodes, x);
  int node_ids_fi2[n_pix] = spline_findpos(fi2_nodes, x);
  int node_ids_logwidth[n_pix] = spline_findpos(logwidth_nodes, x);
  int node_ids_bg[n_pix] = spline_findpos(bg_nodes, x);
  int node_ids_bgsl[n_pix] = spline_findpos(bgsl_nodes, x);
  int node_ids_bgsl2[n_pix] = spline_findpos(bgsl2_nodes, x);
  int node_ids_pm[n_pix] = spline_findpos(pm_nodes, x);
  int node_ids_pmlogwidth[n_pix] = spline_findpos(pmlogwidth_nodes, x);

  matrix[n_pix,4] mat_int = spline_getmat(x, int_nodes,node_ids_int);
  matrix[n_pix,4] mat_fi2 = spline_getmat(x, fi2_nodes, node_ids_fi2);
  matrix[n_pix,4] mat_logwidth = spline_getmat(x, logwidth_nodes, node_ids_logwidth);
  matrix[n_pix,4] mat_bg  = spline_getmat(x, bg_nodes, node_ids_bg);
  matrix[n_pix,4] mat_bgsl = spline_getmat(x, bgsl_nodes, node_ids_bgsl);
  matrix[n_pix,4] mat_bgsl2 = spline_getmat(x, bgsl2_nodes, node_ids_bgsl2);
  matrix[n_pix,4] mat_pm = spline_getmat(x, pm_nodes, node_ids_pm);
  matrix[n_pix,4] mat_pmlogwidth = spline_getmat(x, pmlogwidth_nodes, node_ids_pmlogwidth);

}

parameters
{
  // The parameters of the model are the values at the 
  // spline nodes for the log(intensity) of the stream, 
  // log(width), log(bg), bg slope, bg quadratic slope, 
  // and spatial track
  vector[n_int_nodes] ints;// the faintest is exp(-20)
  vector[n_logwidth_nodes] logwidths;

  vector[n_bg_nodes] bgs;
  vector[n_bgsl_nodes] bgsls;
  vector[n_bgsl2_nodes] bgsl2s;

  vector[n_bg_nodes] bgpm1s;
  vector[n_bgsl_nodes] bgslpm1s;
  vector[n_bgsl2_nodes] bgsl2pm1s;

  vector[n_bg_nodes] bgpm2s;
  vector[n_bgsl_nodes] bgslpm2s;
  vector[n_bgsl2_nodes] bgsl2pm2s;

  vector[n_fi2_nodes] fi2s;
  
  vector[n_pm_nodes] pm1s;
  vector[n_pm_nodes] pm2s;
  vector[n_pmlogwidth_nodes] pmlogwidths;
}

transformed parameters
{
}

model
{


  // 4 locations for var names
  // The evaluation of the spline coefficients 
  // given the values at the spline nodes
  vector[n_int_nodes] coeffs_int= spline_getcoeffs(int_nodes, ints);

  //track 
  vector[n_fi2_nodes] coeffs_fi2= spline_getcoeffs( fi2_nodes, fi2s);

  // width 
  vector[n_logwidth_nodes] coeffs_logwidth= spline_getcoeffs( logwidth_nodes, logwidths);

  // background slopes zero-order in pm
  vector[n_bg_nodes] coeffs_bg = spline_getcoeffs( bg_nodes, bgs);
  vector[n_bgsl_nodes] coeffs_bgsl = spline_getcoeffs( bgsl_nodes, bgsls);
  vector[n_bgsl2_nodes] coeffs_bgsl2 = spline_getcoeffs( bgsl2_nodes, bgsl2s);

  // background slopes 1st order in pm1
  vector[n_bg_nodes] coeffs_bgpm1 = spline_getcoeffs( bg_nodes, bgpm1s);
  vector[n_bgsl_nodes] coeffs_bgslpm1 = spline_getcoeffs( bgsl_nodes, bgslpm1s);
  vector[n_bgsl2_nodes] coeffs_bgsl2pm1 = spline_getcoeffs( bgsl2_nodes, bgsl2pm1s);

  // background slopes 1st order in pm2
  vector[n_bg_nodes] coeffs_bgpm2 = spline_getcoeffs( bg_nodes, bgpm2s);
  vector[n_bgsl_nodes] coeffs_bgslpm2 = spline_getcoeffs( bgsl_nodes, bgslpm2s);
  vector[n_bgsl2_nodes] coeffs_bgsl2pm2 = spline_getcoeffs( bgsl2_nodes, bgsl2pm2s);

  vector[n_pm_nodes] coeffs_pm1 = spline_getcoeffs( pm_nodes, pm1s);
  vector[n_pm_nodes] coeffs_pm2 = spline_getcoeffs( pm_nodes, pm2s);
  vector[n_pmlogwidth_nodes] coeffs_pmlogwidth = spline_getcoeffs( pmlogwidth_nodes, pmlogwidths);
  

  //target+= partial_lpdf_sum(hh,			1,
  target += reduce_sum(partial_lpdf_sum, hh, grainsize,

		       int_nodes,
		       ints,
		       coeffs_int,
		       mat_int,
		       node_ids_int,

		       fi2_nodes,
		       fi2s,
		       coeffs_fi2,
		       mat_fi2,
		       node_ids_fi2,
		       
		       logwidth_nodes,
		       logwidths,
		       coeffs_logwidth,
		       mat_logwidth,
		       node_ids_logwidth,
		       
		       bg_nodes,
		       bgs,
		       coeffs_bg,
		       mat_bg,
		       node_ids_bg,
		       
		       bgsl_nodes,
		       bgsls,
		       coeffs_bgsl,
		       mat_bgsl,
		       node_ids_bgsl,
		       
		       bgsl2_nodes,
		       bgsl2s,
		       coeffs_bgsl2,
		       mat_bgsl2,
		       node_ids_bgsl2,
		       
		       bgpm1s,
		       coeffs_bgpm1,
		       bgslpm1s,
		       coeffs_bgslpm1,

		       bgsl2pm1s,
		       coeffs_bgsl2pm1,
		       bgpm2s,
		       coeffs_bgpm2,
		       
		       bgslpm2s,
		       coeffs_bgslpm2,
		       bgsl2pm2s,
		       coeffs_bgsl2pm2,
		       
		       pm_nodes,
		       pm1s,
		       coeffs_pm1,
		       mat_pm,
		       node_ids_pm,
		       
		       pm2s,
		       coeffs_pm2,
		       
		       pmlogwidth_nodes,
		       pmlogwidths,
		       coeffs_pmlogwidth,
		       mat_pmlogwidth,
		       node_ids_pmlogwidth,
		       
		       x,
		       y,
		       dpm1,
		       dpm2
		       );

  
  target+= normal_lpdf(ints | 0, 20);
  target+= normal_lpdf(logwidths| log(0.8), .3);  
  target+= normal_lpdf(bgs| 0, 100);
  target+= normal_lpdf(bgsls | 0, 100);
  target+= normal_lpdf(bgsl2s| 0, 100);
  target+= normal_lpdf(bgpm1s| 0, 100);
  target+= normal_lpdf(bgslpm1s| 0, 100);
  target+= normal_lpdf(bgsl2pm1s| 0, 100);
  target+= normal_lpdf(bgpm2s| 0, 100);
  target+= normal_lpdf(bgslpm2s| 0, 100);
  target+= normal_lpdf(bgsl2pm2s| 0, 100);
  target+= normal_lpdf(fi2s| 0, 3);
  target+= normal_lpdf(pm1s| 0, 0.5);
  target+= normal_lpdf(pm2s| 0, 0.5);
  target+= normal_lpdf(pmlogwidths| log(.1), .5);

}
