% Tue Jun 19 15:58:39 MSK 2012
% Karl Kästner, Berlin

function fem_plot_1d(name, printflag)
	% load the data
	s = load(name);

	% convergence n vs err_est vs err_true
	err_true = abs(s.E - s.E_true*ones(1,length(s.N(:,2))));
	nErr_true = sqrt(sum(err_true.^2,1));
	figure(1);
	loglog(s.N(:,1),[s.err_est nErr_true']);
	% convergence n vs err_est vs h
	figure(2);
	loglog(s.N(:,1),[s.err_est s.h_side_min']);
	% run time TODO
	% efficency TODO
	% squared eigenvectors
	figure(5);
	[P_ pdx] = sort(s.P);
	for idx=1:s.k
		subplot(ceil(sqrt(s.k)),ceil(sqrt(s.k)),idx);
		plot(P_, s.v(pdx,idx).^2);
	end
end

