% load meshes of different order and get maximum of the histrogram

function func()

F = {'../dat-new/fem-2d-1338937518.mat', '../dat-new/fem-2d-1338937622.mat', '../dat-new/fem-2d-1338937673.mat'};
F = fliplr(F);

clf
%x=linspace(0.0,1.0,1000);
x = logspace(-10,0,1000);

for idx=1:length(F)
		% load the adaptive mesh
		s=load(F{idx});
		%check_v(s)
		[R D] = regularity_2d(s.P, s.T);
		D     = sqrt(sum(D.^2,2));
		err = 2*(2 + 1./D).*R.^(s.order).*exp(-2*D);
		C(1,idx) = max(err)
		r(1,idx) = min(R)
		p(1,idx) = length(s.P)
		%subplot(4,length(F),idx)
		%H = histc(err,x);
		%semilogy(x,H);
		%loglog(x,H);
		subplot(2,length(F),0*length(F)+idx)
		display_2d(s.P,s.T,s.Bc,0);
		axis equal, axis tight
		% reproduce the grid
		tic()
		[P T Bc] = fem_2d_heuristic_mesh(160*[1 1], 160*[0.5 0.5], s.abstol, s.order);
		t(idx) = toc()
		[R D]    = regularity_2d(P, T);
		D        = sqrt(sum(D.^2,2));
		err = 2*(2 + 1./D).*R.^(s.order).*exp(-2*D);
		C(2,idx) = max(err)
		r(2,idx) = min(R)
		p(2,idx) = length(P)
		subplot(2,length(F),1*length(F)+idx)
		display_2d(P,T,Bc,0);
		axis equal, axis tight
end % for

end % function

function check_v(s)
		n=1;
		v = s.v(:,n);
		v = -v(:,n);
		D_ = sqrt(s.P(:,1).^2 + s.P(:,2).^2);
		v_ = 0.5*pi*exp(-2*D_);
		[norm(v) norm(v_ - v)/norm(v)]
		subplot(2,2,1);
		display_2d(s.P,s.T,s.Bc,0, [v(s.T(:,1)) v(s.T(:,2)) v(s.T(:,3))],[],'edgecolor','none');
		axis(10*[-1 1 -1 1]);
		subplot(2,2,2);
		display_2d(s.P,s.T,s.Bc,0, [v_(s.T(:,1),1) v_(s.T(:,2),1) v_(s.T(:,3),1)], [], 'edgecolor','none');
		axis(10*[-1 1 -1 1]);
		subplot(2,2,3);
		v_ = v_ - s.v(:,1);
		display_2d(s.P,s.T,s.Bc,0, [v_(s.T(:,1),1) v_(s.T(:,2),1) v_(s.T(:,3),1)], [], 'edgecolor','none');
		axis(10*[-1 1 -1 1]);
end

