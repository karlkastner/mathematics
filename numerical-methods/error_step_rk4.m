% Tue 12 Mar 16:05:29 CET 2024
function [e,dt0,out] = error_step_rk4(t,z,erel0,final)
	if (nargin()<4 || ~final)
		% compute error for entire time series
		T  = t(end)-t(1);
		nt = length(t);
		dt = T/(nt-1);
		% TODO variable step
	%	D3 = derivative_matrix_k_1d(nt,T,3);
	%	d3z_dt3 = D3*double(z);
	else
		% expects last value in tt(1) and z(:,1)
		A        = vander_1d(t-t(1),5);
		v        = vanderd_1d(0,5,5);
		%d3z_dt3    = v*(z*inv(A)')';
		d5z_dt5    = (z*(A'\v'));
		dt = t(1)-t(2);
		%e        = -0.5*d2z_dt2*dt*dt;
		%e        = e_div_dt2*dt^2;
	end
	e    = (-dt*dt*dt*dt*dt/120)*d5z_dt5;
	rmsz = rms(z(:,1));
	% einf = max(abs(e),[],2);
	erms = rms(e);
	erel = erms./rmsz;
	if (nargin()>2)
		dt0 = dt*(erel0./abs(erel)).^(1/5);
	else
		dt0 = [];
	end
	if (nargout()>2)
		%out.zrms = zrms;
		out.rmsz = rmsz;
		out.erms = erms;
		out.erel = erel;
		out.dt0  = dt0;
	end
end

