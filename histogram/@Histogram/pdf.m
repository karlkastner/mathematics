% 2016-03-04 12:09:23.320601387 +0100
% Karl Kastner, Berlin

function pdf = pdf(hist,x)
	if (nargin()<2)
		x = linspace(hist.edge(1),hist.edge(end));
	end
	pdf = interp1(hist.centre.x,hist.centre.pdf,x,'linear','extrap');
%	pdf = interp1(hist.centre.x,hist.centre.pdf,x,'spline','extrap');
	%pdf(pdf < 0) = 0;
	pdf(x < hist.edge.x(1) | x>hist.edge.x(end)) = 0;
	if (nargout()<1)
		plot(x,pdf);
	end
end

