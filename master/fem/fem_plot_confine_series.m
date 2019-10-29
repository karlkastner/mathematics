% Sat Jun 23 18:12:28 MSK 2012
% Karl Kästner, Berlin


function fem_plot_confine_series(folder,mode, printflag)
	bw = 0;
	name_c=regexp(ls(folder,'-1'),'\n','split')
	folder_c = regexprep(folder, '/$', '');
	folder_c = regexprep(folder_c, '.*/', '');
	ms = 12;
	lw = 2;
	%folder_c=regexp(folder,'/','split')
	if (nargin() < 2 || isempty(mode)) mode = 1; end

	idx_ = 1; 
	for idx=1:length(name_c)
		% load data
		try
		s = load([folder '/' name_c{idx}]);
		if (~isfield(s,'convflag')) s.convflag=0; end
		s=orderfields(s);
		S(idx_) = s;
		% extract some values
		s = S(idx_);
		x0(idx_) = s.x0(1);
		L0(idx_) = s.L0(1);
		Tr(idx_) = sum(sum(s.Tr));
		N(idx_) = s.N(end,1);
		E(:,idx_) = s.E(:,end);
		idx_ = idx_+1;
		catch exception
			exception
			'failed to load'
			name_c{idx}
		end
	end

	% sort for increasing domain size
	figure(1); clf
	switch (mode)
		case{1}
		[L0 idx] = sort(L0);
		E  = E(:,idx);
		Tr = Tr(:,idx);
		N  = N(:,idx);
		S  = S(:,idx);
		l  = length(L0);

		% plot the eigenvalues
		%subplot(1,2,2);
		semilogx(L0,E','.-k');
		ylim([-2 10]);
		xlim([5e-1 5e1]);
%		axes('pos',[0.5 0.6 0.4 0.3]);
%		semilogx(L0,E','.-k');
%		ylim([-2 15]);
%		xlim([5e-1 1e2]);
%		semilogx(L0,E','.-');
%		loglog(L0,E'+3,'.-');
		%plot(L0,E','.-');
		%plot(sqrt(L0),sign(E).*sqrt(abs(E)),'.-');
		%semilogx(L0,sign(E).*log(abs(E)),'.-');
		%semilogx(L0,E,'.-');
		xlabel('L [a_0]')
		ylabel('\lambda')
		grid on
		%xlim([1 32])
		% plot number of gridpoints
		figure(2);
		subplot(1.5,1.5,1);
		loglog(L0,N,'k.-','linewidth',lw,'markersize',ms);
		xlabel('L [a_0]');
		ylabel('n : number of gridpoints')
		grid on
		%set(gca,'minorgrid','none');
		% plot run time
		figure(3);
		subplot(1.5,1.5,1);
		loglog(L0,Tr,'k.-','linewidth',lw,'markersize',ms);
		xlabel('L [a_0]');
		ylabel('run time [s]');
		grid on;
		%set(gca,'minorgrid','none');

		case{2}
		[x0 idx] = sort(x0);
		E = E(:,idx);
		Tr = Tr(:,idx);
		N  = N(:,idx);
		S = S(:,idx);
		l = length(x0);

		figure(1);
%		loglog(x0,E'+3,'.-');
		semilogx(x0,E,'.-k')
		xlabel('x_0 [a_0]');
		ylabel('\lambda')
		ylim([-0.3 0.1])
		grid on
		

		% plot number of gridpoints
		figure(2);
		loglog(x0,N,'.-');
		xlabel('x_0 [a_0]');
		ylabel('number of gridpoints')
		grid on

		% plot run time
		figure(3);
		loglog(x0,Tr,'.-');
		xlabel('x_0 [a_0]');
		ylabel('run time [s]')
		grid on
	end

	if (nargin > 2 && printflag)
		figure(1)
		preparePrint();
		name_ = ['../img/' folder_c '-eigenvalues.eps'];
		name_
		print('-deps', name_);
		system(['epstopdf ' name_]);

		figure(2)
		preparePrint();
		name_ = ['../img/' folder_c '-convergence.eps'];
		print('-deps', name_);
		system(['epstopdf ' name_]);

		figure(3)
		preparePrint();
		name_ = ['../img/' folder_c '-run-time.eps'];
		print('-deps', name_);
		system(['epstopdf ' name_]);
	end

	% generate a table
	if (1 == mode)
		for idx=1:length(L0)
			if ismember( single(L0(idx)), single(2.^(-5:0.5:8)))
			fprintf('%7.3f ',L0(idx));
			fprintf(1,'%+1.3e ',E(:,idx)');
			fprintf(1,'\n');
			end
		end
	end
	if (2 == mode)
		for idx=1:length(x0)
			if ismember( single(x0(idx)), single(2.^(-5:0.5:8)))
				fprintf('%7.3f ',x0(idx));
				fprintf(1,'%+1.3e ',E(:,idx)');
				fprintf(1,'\n');
			end
		end
	end

	return;
	ss = 0;
	F_old = [];
	figure(2);
	for idx=1:1 %size(E,1)
%		figure(idx+1); clf
%		subplot(ceil(sqrt(length(E(:,1)))),ceil(sqrt(length(E(:,1)))),idx);
		for jdx=1:l
		figure(jdx+1); clf
			s = S(jdx);
			mesh = Mesh(s.P,s.T,s.Bc);
			lc = ceil(sqrt(l));
			lr = ceil(l/lc);
%			subplot(lr,lc,jdx);
			%subplot(round(sqrt(size(E,1))),round(sqrt(size(E,1))),jdx);
			cla();
%			pause
			subplot(2,2,1)
			display_2d(mesh, 0, [s.v(s.T(:,1),idx).^2 s.v(s.T(:,2),idx).^2 s.v(s.T(:,3),idx).^2], [], 'EdgeColor', 'k');
			axis equal; axis tight

			[area l_boundary h_side s_angle C] = regularity_2d(s.P,s.T,s.Bc);
			Nm = double(FEM.elementNeighbours(s.P,s.T,s.Bc));
			mesh = Mesh(s.P,s.T,s.Bc);
			[M v_err err_max nH dV] = mark_2d_10(mesh, Nm, s.v(:,1), area, h_side, s_angle, C);
			d1V = fem_2d_d1V(s.T,s.P,s.v(:,1),0);
			%v = v_err; %s.v_err;
			%v = log(v);
			%d1V = (sqrt(sum(d1V.^2,2)));
			%d1V = (sqrt(d1V(:,1).^2 + d1V(:,2).^2));
			%d1V = d1V.^2;
			%v = d1V;
			s.opt.backward
			v = v_err;
			subplot(2,2,2)
			display_2d(mesh, 0, v, [], 'EdgeColor', 'k');
			%display_2d(mesh, 0, [v(s.T(:,1),idx) v(s.T(:,2),idx) v(s.T(:,3),idx)], [], 'EdgeColor', 'k');
			axis equal; axis tight
			colorbar
			%v = s.v(:,1);

			%v = dV;
			%v = d1V.*prod(h_side(:,1),2).^(1/3);
			v = prod(h_side(:,1),2).^(1/3);
			subplot(2,2,3)
			display_2d(mesh, 0, v(:,1), [], 'EdgeColor', 'k');
			%display_2d(mesh, 0, [s.v(s.T(:,1),idx).^2 s.v(s.T(:,2),idx).^2 s.v(s.T(:,3),idx).^2], [], 'EdgeColor', 'k');
			colorbar
			axis equal; axis tight

			v = s.v_err;
			subplot(2,2,4)
			display_2d(mesh, 0, v, [], 'EdgeColor', 'k');
			colorbar
			axis equal; axis tight

			%subplot(2,2,3)
			%display_2d(mesh, 0, [s.v(s.T(:,1),idx).^2 s.v(s.T(:,2),idx).^2 s.v(s.T(:,3),idx).^2], [], 'EdgeColor', 'k');
			%axis equal; axis tight
%			s
			%v = log(s.v_err);


			%v = s.v;
%			size(dV)
%			dV = fem_2d_d3V(s.T, s.P, s.v(:,1))
%			v = dV(:,1);
%			subplot(1,2,1)
			%display_2d(mesh, 0, [v(s.T(:,1),idx).^2 v(s.T(:,2),idx).^2 v(s.T(:,3),idx).^2], [], 'EdgeColor', 'k');
%			display_2d(mesh, 0, v, [], 'EdgeColor', 'k');
%			axis equal; axis tight
%			colorbar
%			subplot(1,2,2)
%			v = s.v(:,1);
%			display_2d(mesh, 0, [v(s.T(:,1),idx).^2 v(s.T(:,2),idx).^2 v(s.T(:,3),idx).^2], [], 'EdgeColor', 'k');
%			axis equal; axis tight
%			colorbar
%			pause

			if (bw)
				colormap gray
				F = getframe();
				n = 128;
				data = imresize(F.cdata(:,:,1),[n,n]);
				level = graythresh(data);
				data_bw = im2bw(data,level);
				imagesc(data_bw)
			end

			%title([num2str(ss(jdx))])
			title(s.E(idx));
			axis equal, axis tight
		end
		x0(jdx)
	end		
end % fem_plot_confine_series()

