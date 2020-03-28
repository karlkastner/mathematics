% Tue May 15 16:56:47 MSK 2012
% Karl Kästner, Berlin

function fem_plot_2d(name, printflag)

	if (nargin <2)
		printflag = 0;
	end

	lw=2;
	ms=12;
	
        dirname  = regexp(name, '.*/', 'match')
        basename = regexprep(name, '.*/', '')

	% load the data
	s = load(name);
	%if (~isfield(s,'mesh'))
	s.mesh = Mesh_2d(s.P,s.T,s.Bc);
	if (isfield(s,'h_side_min'))
		s.h_side = s.h_side_min;
	end
	if (isfield(s,'h_min'))
		s.h_side = s.h_min;
	end
	if (isfield(s,'s_angle_min'))
		s.degen = s.s_angle_min;
	end

	s.K = s.K(1:size(s.E,2));
	s.N = s.N(1:size(s.E,2),:);
	s.Tr = s.Tr(1:size(s.E,2),:);

	%disp(sprintf('%s order %d L0 %s x0 %s k %d abstol %e', name{idx}, order, mat2str(L0), mat2str(x0), k, abstol ));
	% plot convergence
	figure(1); clf; subplot(1.5,1.5,1);
	if (isempty(s.E_true))
		s.E_true = s.E(:,end);
	end
	Err = abs(s.E - s.E_true*ones(1,size(s.E,2)));
	% Err = abs(s.E - s.E_true*ones(1,length(s.N(:,2))));
	 %./(s.E_true*ones(1,length(s.N(:,2)))));
	nErr = sqrt(sum(Err.^2,1));
	%[ nErr' s.err_est h_side [0 diff(E(1,:))]' ]'
	switch (printflag)
		case {1}
			loglog(s.N(:,1), [ s.err_est nErr' ],'.-');
		case {2}
			loglog(s.N(:,1), [ s.err_est nErr' s.h_side [0 abs(diff(E(1,:)))]'*0.5^(1) MM],'.-');
		otherwise
			size(s.N,1)
			size(nErr)
			size(s.err_est)
			loglog(s.N(:,1), [ s.err_est nErr' ],'.-');
	end
	grid on; set(gca,'minorgrid','none');
	title('Absolut Discretisation Error');
	xlabel('total number of grid points : n');
	legend('location','southwest','estimated ||v-v_*||_{L\infty}', 'true |\lambda-\lambda_*|/|\lambda_*|');

	% Mesh Quality
	figure(8); clf; subplot(1.5,1.5,1);
	loglog(s.N(1:length(nErr),1), [nErr' s.h_side' s.degen'],'.-','linewidth',lw,'markersize',ms);
	grid on;
	set(gca,'minorgrid','none');
	xlim([1e3 1e5])
	xlabel('total number of grid points : n');
	legend('location','southwest',{'estimated $||v {-} v_*||_{L_\infty}$', '$h_{min}$', '$\sin(\alpha_{min})$'},'interpreter','latex');

	% plot convergence only of the latest eigenvalue
	%plotyy()
	figure(2); clf; subplot(1.5,1.5,1);
	%if (exist('K','var'))
	if (isfield(s,'K'))
		loglog( s.N(:,1), [s.err_est Err(sub2ind(size(Err),s.K',1:size(Err,2)))'],'.-');
%		loglog( N(:,1), [[s.err_est(1); s.err_est(1:end-1)] Err(sub2ind(size(Err),K',1:size(Err,2)))'],'.-');
%		loglog( N(:,1), [s.err_est Err(sub2ind(size(Err),min(K(end),K'+1),[2:size(Err,2) size(Err,2)])')],'.-');
		grid on; set(gca,'minorgrid','none');
		title('Relative Discretisation Error of Latest Eigenvalue');
		xlabel('total number of grid points : n');
		legend('location','southwest','estimated ||v-v_*||_{L\infty}', 'true |\lambda-\lambda_*|/|\lambda_*|');	
	end

	% plot run time
	figure(3); clf; subplot(1.5,1.5,1);
%	if (~(exist('s.Tr','var')))
%		s.Tr = Ti;
%	end
	Tr = s.Tr; Tr(end,end) = NaN;
	if (size(Tr,2) > 3) Tr=Tr(:,2:4); end
	loglog(s.N(:,1),Tr,'.-','linewidth',lw,'markersize',ms);
	%loglog(s.N(:,1),[s.Tr sum(s.Tr,2) cumsum(sum(s.Tr,2))],'.-');
	grid on; set(gca,'minorgrid','none');
	%title('Run Time');
	xlabel('total number of grid points : n');
	ylabel('time [s]')
	if (size(Tr,2) < 4)
		legend('location','northwest','assembly','eigenvalue','refinement'); %,'last iteration total','accumulated')
	else
		legend('location','northwest','assembly','eigenvalue','refinement'); %,'last iteration total','accumulated')
		%legend('location','northwest','mesh generation', 'assembly','eigensolver','refinement'); %,'last iteration total','accumulated')
	end
	xlim([1e3 1e5])

	% plot efficiency
	figure(4); clf; subplot(1.5,1.5,1);
	loglog(cumsum(sum(s.Tr,2)), nErr','.-');
	title('Efficency');
	grid on; set(gca,'minorgrid','none');
	xlabel('time [s]')
	ylabel('|\lambda - \lambda_*|/|\lambda_*|');

	% plot the eigenvector(s)
	figure(5); clf();
	n = size(s.v,2);
	rows = round(sqrt(n));
	cols = ceil(n/rows);
	for idx=1:n
		subplot(rows,cols,idx);
		display_2d(s.mesh, 0, [s.v(s.T(:,1),idx).^2 s.v(s.T(:,2),idx).^2 s.v(s.T(:,3),idx).^2], [], 'EdgeColor', 'none');
		title(['\lambda_{' num2str(idx) '} = ' num2str(s.E(idx,end))]);
		axis equal; axis tight;
	end
	
	% plot the estimated error (only precomputed for final vector)
	figure(6); clf
	display_2d(s.mesh,0,log10(s.v_err),[],'EdgeColor','none');
	title('Interpolation Error');
	axis equal; axis tight; colorbar;
	
	% plot the mesh
	figure(7); clf
	display_2d(s.mesh,0,[],[],'facecolor','w');
	axis equal; axis tight

	% plot the radial wave functions for each eigenvector
	% this is a rudimentary implementation
	% TODO : implement by grid sampling and not by summing triangles
	figure(9); clf
	n_bin = 32;
	% TODO, what if not centred
	r_max = max(s.L0)*sqrt(0.5);
	for vdx=1:size(s.v,2)
		[R Vr(:,vdx)] = fem_get_2d_radial(s.P, s.T, s.v(:,vdx), n_bin, r_max);
	end % for vdx
	%plot(R, Vr);
	for idx=1:n
		subplot(rows,cols,idx);
		plot(R, Vr(:,idx));
		axis([0 r_max 0 1]);
	end

	% plot degenerated triangles	
%	figure(1); clf
%	R = regularity_2d(s.P, s.T);
%	display_2d(s.P,s.T,s.Bc,0,log10(R));
%	axis equal; axis tight
%	colorbar;
%	title('Triangle Regularity and Degeneracy')

	s.opt.poly
	s.convflag
	max(s.N(:,1))
	fprintf(1,'L0: %f %f\n', s.L0);
	fprintf(1,'x0: %f %f\n', s.x0);
	fprintf(1,'reltol: %1.0e\n', s.opt.reltol);
	fprintf(1,'abstol: %1.0e\n', s.opt.abstol);
	fprintf(1,['%2d %+1.3e\n'],[ (1:size(s.E,1))' s.E(:,end)]');
	
	if ( 0 ~= printflag )
		name_ = [ '../img/' basename '-convergence.eps'];  figure(1); preparePrint(); print('-depsc', name_);
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-run-time.eps'];     figure(3); preparePrint(); print('-depsc', name_);
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-efficency.eps'];    figure(4); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-eigenvector.eps'];  figure(5); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-errorvector.eps'];  figure(6); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-mesh.eps'];         figure(7); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-mesh-quality.eps']; figure(8); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
		name_ = [ '../img/' basename '-radial-wave-functions.eps']; figure(9); preparePrint(); print('-depsc', name_); 
		system(['LD_LIBRARY_PATH= epstopdf ' name_]); system(['rm ' name_]);
	end
end % fem_plot_2d()

