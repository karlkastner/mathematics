% Mon Jul 30 17:26:59 MSK 2012
% Karl Kästner, Berlin

function plot_fem_3d_series(name_c, folder, vecflag, printflag)

	if (nargin() < 2)
		printflag = 0;
	end

	if (isempty(name_c))
		name_c=ls([folder '/*.mat'],'-1');
%		name_c=ls(folder,'-1',' *.mat');
		name_c=name_c(1:end-1); % strip newline
		name_c=regexp(name_c,'([^ \t\n]*)','match');
	end

	% clear the figures
	for idx=1:7
		figure(idx); clf();
	end
	
	lw=2;
	ms=12;
	colour = {'b','g','r','k','m','y', 'c', 'r','g'};
	n_row = 2;
	n_col = ceil(length(name_c)/n_row);
	
	for idx=1:length(name_c)
		name_c{idx} = regexprep(name_c{idx},'.*/','')
		
		% load data
%		try
		if (~isempty(folder)) folder = [folder '/']; end
		s = load([folder name_c{idx}]);
		s = load([folder name_c{idx}]);
		mesh = Mesh_3d(s.P,s.T,s.Bc);
		
		% convergence
		figure(1); %if (0 == printflag) subplot(1,4,1); else
		subplot(1.5, 1.5, 1); 
		Err = (s.E - s.E_true*ones(1,length(s.E(:,2))))./(s.E_true*ones(1,length(s.E(:,2))));
		nErr = sqrt(sum(Err.^2,1));
		disp(s.opt.poly);
		%nErr'./err_est
		loglog(s.N(1:size(Err,2),1),nErr,'.-','color',colour{idx},'linewidth',lw,'markersize',ms);
		hold on;
		grid on; set(gca,'minorgrid','none');
		%title('Convergence');
		xlabel('total number of grid points : n');
		ylabel('|\lambda - \lambda_*|/|\lambda_*|');
		xlim([1e2 1e5])
		
		figure(6)
		loglog(s.N(1:size(s.h_min,2),1),s.h_min,'.-','color',colour{idx}); hold on;
		loglog(s.N(1:size(s.degen,2),1),s.degen,'o-','color',colour{idx}); hold on;
		loglog(s.N(1:size(nErr,2),1),nErr,'.--','color',colour{idx}); hold on;
		grid on; set(gca,'minorgrid','none');
		xlabel('total number of grid points : n');
		ylabel('|\lambda - \lambda_*|/|\lambda_*|');
		xlim([1e2 1e5])

		figure(7)
%		log(nErr'./err_est)
%		C_ = [1e-0, 1e-2, 1e-4, 1e-6, 1e-8];
%		err_est = err_est/C_(opt.poly);
%		C__ = [-1.0736, -4.5535, -8.8459, -12.4002];
%		nErr = nErr.*h_min.^(-2*opt.poly + opt.poly + 1 );
%		mean(log10(nErr'./err_est))
%		err_est = err_est*10^C__(opt.poly);
%		mean(log10(nErr'./err_est))
		loglog(s.N(1:length(s.err_est),1),s.err_est,'o-','color',colour{idx}); hold on;
		loglog(s.N(1:size(nErr,2),1),nErr,'.--','color',colour{idx}); hold on;
		grid on; set(gca,'minorgrid','none');
		xlabel('total number of grid points : n');
		ylabel('|\lambda - \lambda_*|/|\lambda_*|');
		xlim([1e2 1e5])

		%subplot(1,4,1);
		%loglog(s.N(:,1),err_est,'.--','color',colour{idx}); hold on;
		%grid on; set(gca,'minorgrid','none');
		%title('Convergence'); %Relative Discretisation Error');
		%xlabel('total number of grid points : n');
		%ylabel('|\lambda - \lambda_*|/|\lambda_*|');
		%figure(1)
		
		% run time
		figure(2);
		%if (0 == printflag) subplot(1,4,3); else figure(2); subplot(1.5, 1.5, 1); end
		loglog(s.N(1:size(s.Tr,1),1), cumsum(sum(s.Tr,2)), '.-', 'color', colour{idx},'linewidth',lw,'markersize',ms);
		hold on;
		grid on; set(gca,'minorgrid','none');
		title('Run time');
		xlabel('total number of grid points : n');
		ylabel('time [s] (accumulated)');

		% Efficency
		figure(3);
		%if (0 == printflag) subplot(1,4,4); else figure(3);
		subplot(1.5, 1.5, 1);
		loglog(cumsum(sum(s.Tr(1:length(nErr),:),2)),nErr,'.-','color',colour{idx},'linewidth',lw,'markersize',ms);
		hold on
		grid on; set(gca,'minorgrid','none');
		%title('Efficency');
		xlabel('time [s] (accumulated)');
		ylabel('|\lambda - \lambda_*|/|\lambda_*|');
		%xlim([0.5 200])


		if (nargin > 2 && vecflag)

		% error vector
		figure(4);
		subplot(n_row,n_col,idx);
		display_2d(mesh,0,log10(v_err),[],'EdgeColor','none'); hold on
		caxis([-5.5 -3]); colorbar; axis equal; axis tight;
		title(['Basis ' num2str(opt.poly)']);

		% ground states vectors
		figure(5);
		subplot(n_row,n_col,idx);
		display_2d(mesh, 0, [s.v(s.T(:,1)).^2 s.v(s.T(:,2)).^2 s.v(s.T(:,3)).^2], [], 'EdgeColor','none'); hold on
		colorbar; axis equal; axis tight;
		title(['Basis ' num2str(s.opt.poly)'])

		end
		
		L(idx,1) = s.opt.poly;
%		L_(2*idx-1,1) = order;
		KK(idx,1) = s.k;
%		catch err
%			err
%		end
	end % for idx

	% plot legends
	figure(1); %if (0 == printflag) subplot(1,4,1); end;
	legend('location', 'southwest', [num2str(L, 'p=%d') num2str(KK,' k=%d')]);
	figure(2); % if (0 == printflag) subplot(1,4,2); else figure(2); end
	legend('location', 'northwest', [num2str(L, 'p=%d') num2str(KK,' k=%d')]);
	figure(3); %if (0 == printflag) subplot(1,4,3); else figure(3); end
	legend('location', 'southwest', [num2str(L, 'p=%d') num2str(KK,' k=%d')]);
	figure(6)
	legend('location', 'southwest', [num2str(L, 'p=%d') num2str(KK,' k=%d')]);
	%figure(4); %if (0 == printflag) subplot(1,4,3); else figure(3); end
	%legend('location', 'southwest', [num2str(L, 'Poly=%d') num2str(KK,' k=%d')]);

	if (nargin()>3 && printflag)
		name_ = [ '../img/' name_c{1} '-series-convergence.eps'], figure(1); preparePrint();  print('-depsc', name_); system(['epstopdf ' name_]);
		name_ = [ '../img/' name_c{1} '-series-run-time.eps'];    figure(2); preparePrint();  print('-depsc', name_); system(['epstopdf ' name_]);
		name_ = [ '../img/' name_c{1} '-series-efficency.eps'];   figure(3); preparePrint();  print('-depsc', name_); system(['epstopdf ' name_]);
	end
end % function plot_fem_2d_series()

