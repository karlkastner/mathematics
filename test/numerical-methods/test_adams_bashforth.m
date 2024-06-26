% test adams bashforth
	%f' = x

	df_dt = @(t,f) f;
	df_dt = @(t,f) (2*pi).*cos(2*pi*t);
%	df_dt = 
	f0 = 0;
%	f0 = [1];

%	% exponential function
%	eps = 
	T = [0,2];

	nt = 100*T(end);
	ks = 1;
	f = [];
	for norder=1:4
		[t,f(:,norder)] = ode_adams_bashforth(df_dt,f0,T,nt,norder,ks);
	end
	%f(:,end+1) = exp(t);
	f_ = sin(2*pi*t); 
	res = f-f_;
%(:,end);
	rmsres = rms(res);

	figure(1);
	clf();
	subplot(2,2,1)
	plot(t,f);
	hold on
	plot(t,f_,'k--')
	legend
	subplot(2,2,2)
	plot(t,res)

	subplot(2,2,3)
	semilogy(rmsres)
	


	

