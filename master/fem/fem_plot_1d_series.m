% Tue Jun 19 17:17:53 MSK 2012
% Karl Kästner, Berlin


function plot_series(folder, prinflag)
	name_c=regexp(ls(folder,'-1'),'\n','split');
 
	% load data sets
	for idx=1:length(name_c)
		try
			S(idx) = load([folder '/' name_c{idx}]);
			s = S(idx);
			L0(idx) = s.L0(1);
			x0(idx) = s.x0(1);
			E(:,idx) = s.E(:,end);
		catch
		end
	end
	figure(2); clf
	%semilogx(L0,E,'.-');
	k = 4;
	%loglog(L0,E(1:k,:)+1.5,'.-');
	semilogx(L0,E(1:k,:),'.-');
	xlabel('L [a_0]');
	ylabel('\lambda+3');
	xlim([1 32]);
	ylim([-20 140]);
	grid on;

	if (nargin() > 1 || && printflag)
		preparePrint();
%		name = ['img/1d-confine-series.
	end
end % plot_series

