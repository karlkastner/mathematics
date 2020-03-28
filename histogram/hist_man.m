% 2015-02-05 19:52:35.156258330 +0100
% Karl Kastner, Berlin

% TODO, smooth by piecewise cubic polynomials
% TODO this is and ad-hoc implementation
%      established methods can be found in the literature under
%      "kernel density estimation" (KDE)
%
% with n_elem = sqrt(n_x)/4
function [q_ pdf_ cdf_] = hist_man(x,n,varargin)
	fdx = isfinite(x);
	m = size(fdx,1);
	if (nargin < 2 || isempty(n))
		n = ceil(sqrt(m));
	else
		n = min(n,m);
	end
	cdf = (1:n)'/(n+1);
	cc  = colormap('lines');
	ih = ishold();
	ci = get(gca, 'ColorOrderIndex');
	if (ih)
		v=ver('matlab');
		if (datenum(v.Date) < datenum('01/01/2014'))
			ci = ci+1;
		end
	end
%	if (~ih)
%	cla();
%	end
	for idx=1:size(x,2)
		q = quantile(x(fdx(:,idx)),cdf);
		[q id] = unique(q);
		q0 = linspace(q(1),q(end),n)';
		pdf = diff(cdf(id))./diff(q);
		pdf(end+1) = pdf(end);
		cdf = interp1(q,cdf(id),q0);
		pdf = interp1(q,pdf,q0);
		q = q0;
		if (0 == nargout())
			if (0 == length(varargin))
				set(gca, 'ColorOrderIndex', ci);
				plot(q,pdf,varargin{:});
				v=ver('matlab');
				if (datenum(v.Date) < datenum('01/01/2014'))
					ci = ci+1;
				end
			else
				plot(q,pdf,varargin{:});
			end
			hold on
		end
		q_(:,idx) = q;
		pdf_(:,idx) = pdf;
		cdf_(:,idx) = cdf;
	end
	if (~ih)
		hold off
	end
end % hist_man

