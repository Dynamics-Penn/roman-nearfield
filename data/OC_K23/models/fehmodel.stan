data
{
  int nsel;
  vector[nsel] feh_sel;
  vector[nsel] fi1_sel;
  vector[nsel] efeh_sel;

  int nbg;
  vector[nbg] feh_bg;
  vector[nbg] efeh_bg;
  vector[nbg] fi1_bg;

  real grad_prior_width;
}

parameters
{
  real feh_mean_stream;
  real feh_mean_bg;

  real grad_stream;
  real grad_bg;

  real<lower=0> feh_sig_stream;
  real<lower=0> feh_sig_bg;

  real<lower=0,upper=1> frac_sel;
  real<lower=0,upper=1> frac_bg;

}

model
{
  grad_stream ~ normal(0, grad_prior_width);
  grad_bg ~ normal(0, grad_prior_width);

  feh_mean_stream ~ normal (-2,2);
  feh_mean_bg ~ normal (-2,2);

  feh_sig_bg ~ normal(0,2);
  feh_sig_stream ~ normal(0,2);
  
  target += -  10000 * fdim(frac_bg, frac_sel);
  
  for (i in 1:nsel)
    {
      target += log_mix(frac_sel,
			normal_lpdf(feh_sel[i]
				    | feh_mean_stream +
				    fi1_sel[i] * grad_stream,
				    sqrt(feh_sig_stream^2 + efeh_sel[i]^2)),
			normal_lpdf(feh_sel[i]| feh_mean_bg +
				    fi1_sel[i]  * grad_bg,
				    sqrt(feh_sig_bg^2 + efeh_sel[i]^2)));
			
    }
  for (i in 1:nbg)
    {
      target += log_mix(frac_bg,
			normal_lpdf(feh_bg[i]
				    | feh_mean_stream +
				    fi1_bg[i]  * grad_stream,
				    sqrt(feh_sig_stream^2 + efeh_bg[i]^2)),
			normal_lpdf(feh_bg[i]| feh_mean_bg +
				    fi1_bg[i]  * grad_bg,
				    sqrt(feh_sig_bg^2 + efeh_bg[i]^2)));
			
    }
  
}
