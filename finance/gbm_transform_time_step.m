% Sun 12 Jan 11:11:30 +08 2020
function [r_,sigma_] = gbm_transform_time_step(r,sigma,dt_old,dt_new)
	r_     = (1+r)^(dt_old/dt_new) - 1;
	sigma_ = sigma*sqrt(dt_old/dt_new);
end


