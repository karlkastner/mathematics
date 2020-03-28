% Tue May 22 19:42:02 MSK 2012
% Karl Kästner, Berlin

function fdm_plot(name_c, pflag)
	% do not get screwed up by buggy ATI/AMD drivers in Linux
	opengl neverselect;

	if (nargin() < 2)
		pflag = 0;
	end

	% reset figures
	for idx=1:7
		figure(idx); clf
	end

	colour = {'b','g','r','k','m'};

	% for all input files
	for jdx=1:length(name_c)
	
		load(name_c{jdx});
		X1 = X{1};
		n1 = size(X1,1);
		if (d > 1)
			X2 = X{2};
			n2 = size(X2,1);
		end
		if (d > 2)
			X3 = X{3};
			n3 = size(X3,1);
			I1 = ones(n1,1);
			I2 = ones(n2,1);
			I3 = ones(n3,1);
		end
	
	        if (isempty(E_true))
	                err_true = [];
	        else
		        err_true = abs ( ( E - E_true*ones(1,size(E,2))).^2 ...
					./ (E_true*ones(1,size(E,2))));
			nerr_true = sqrt(sum(err_true,1));
	        end
	
		% Convergence
		figure(1); subplot(1.5,1.5,1);
		loglog(N, [err_est' nerr_true'],'.-'); hold on;
		grid on; set(gca,'minorgrid','none');
		title('Convergence');
		xlabel('total number of grid points : n')
		legend('estimated ||v - v_*||_{L1}', 'true |\lambda - \lambda_*|/|\lambda_*|');
	
		% Convergence of individual eigenvalues
		figure(2); subplot(1.5,1.5,1);
		% TODO - only latest eigenvalue from K(1:end)
		loglog(N, [err_true'],'.-b'); hold on;
		grid on; set(gca,'minorgrid','none');
		xlabel('total number of grid points : n')
		title('Individual Eigenvalues');
	
		% run time
		figure(3); subplot(1.5,1.5,1);
		loglog(N, [cumsum(sum(Tr))'],'.-','color', colour{jdx}); hold on
		grid on; set(gca,'minorgrid','none');
		title('Run Time');
		xlabel('total number of grid points : n')
		ylabel('time [s]')
	
		% efficency
		figure(4); clf;
		loglog(cumsum(sum(Tr)), [err_est' nerr_true'],'.-'); hold on;
		grid on; set(gca,'minorgrid','none');
		title('Effciency');
		xlabel('time [s] (accumulated)')
		legend('estimated ||v - v_*||_{L1}', 'true |\lambda - \lambda_*|/|\lambda_*|');
	
		% eigenvector(s)
		figure(5);
		n = size(v,2);
		rows = round(2*sqrt(n/6));
		cols = ceil(n/rows);
		for idx=1:n
			subplot(rows,cols,idx);
			if (2 == d)
				v_ = reshape(v(:,idx), size(X{2},1)-2, size(X{1},1)-2);
				v_ = padarray(v_,[1 1]);
				[XX YY] = meshgrid(X{1},X{2});
				h=surf(XX,YY,v_.^2); shading(gca,'interp');
				%set(h, 'edgecolor','k');
				view([0 90]);
			else
				v_ = reshape(v(:,idx), n3-2, n2-2, n1-2);
				v_ = padarray(v_.^2,[1 1 1]);
				[XX YY ZZ] = meshgrid(X1,X2,X3);
				v_ = shiftdim(v_,1);
				h = slice(XX, YY, ZZ, v_, 0, 0 ,0);
	%			set(h,'EdgeColor', 'none');
			end
				axis square; axis equal, axis tight
				title(['\lambda_{' num2str(idx) '} = ' num2str(E(idx,end))]);
		end % for idx
	
		% error (only precomputed for last vector)
		figure(6); subplot(1.5,1.5,1);
		if (2 == d)
			v_ = log10(reshape(err, n2-1, n1-1));
			%v_ = padarray(v_,[1 1]);
			X1c = 0.5*(X1(1:end-1) + X1(2:end));
			X2c = 0.5*(X2(1:end-1) + X2(2:end));
			[XX YY] = meshgrid(X1c,X2c);
			h = surf(XX,YY,v_);
		        shading(gca,'interp');
		%	set(h, 'edgecolor','k');
			view([0 90])
		else
			v_ = log10(err);
			size(err)
			% the error is already pre-shifted and dimensionalised
			X1c = 0.5*(X1(1:end-1) + X1(2:end));
			X2c = 0.5*(X2(1:end-1) + X2(2:end));
			X3c = 0.5*(X3(1:end-1) + X3(2:end));
			[XX YY ZZ] = meshgrid(X1c,X2c,X3c);
			h = slice(XX, YY, ZZ, v_, 0, 0, 0);
	%		set(h,'EdgeColor', 'none');
		end
		axis square; axis equal, axis tight
		title('Error Estimate of the Eigenvector')
		colorbar();
	
		% Grid
		figure(7); subplot(1.5,1.5,1);
		if (2 == d)
			% reshape eigenvector and pad boundary values
			v_ = reshape(v(:,1), n2-2, n1-2);
			v_ = padarray(v_,[1 1]);
			[XX YY] = meshgrid(X1,X2);
			h=surf(XX,YY,v_,'facecolor','w');
		%      	shading(gca,'interp');
		%	set(h, 'edgecolor','k');
			view([0 90]);
		else
			v_ = reshape(v(:,1), n3-2, n2-2, n1-2);
			v_ = padarray(v_,[1 1 1]);
			[XX YY ZZ] = meshgrid(X1,X2,X3);
			v_ = shiftdim(v_,1);
			% TODO, white face
			h = slice(XX, YY, ZZ, 0.*v_, 0, 0, 0);
			set(h,'EdgeColor', 'none');
		end
		axis square; axis equal, axis tight
		colorbar();
		title('Eigenvector')
	
		KK(jdx,1) = k;
	end % for jdx

	figure(3);
	legend('location', 'NorthWest', num2str(KK,'k = %d'));

	% print the images to files
	if (pflag)
		name = name_c{1};
		figure(1);
		print([name '-series-convergence.eps'],'-depsc'); 
		figure(3);
		print([name '-series-run-time.eps'],'-depsc'); 
		figure(4);
		print([name '-series-efficency.eps'],'-depsc'); 
		figure(5);
		print([name '-series-eigenvector.eps'],'-depsc'); 
		figure(6);
		print([name '-series-errorvector.eps'],'-depsc'); 
		figure(7);
		print([name '-series-grid.eps'],'-depsc'); 
		system('a=$(ls *.eps); for i in $a; do epstopdf $i; done');
	end
end % function fdm_plot_2d()

