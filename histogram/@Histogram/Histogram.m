% Di 8. Mär 09:24:38 CET 2016
% Karl Kastner, Berlin
classdef Histogram < handle
	properties
		edge
		h
		method = 'constant';
	end % properties
	properties (Constant)
		SHEPPARD = 0;		
		METHOD   = 'constant';		% 2
		%METHOD  = 'midpoint';		% 1 not strictly positive var
		%METHOD  = 'trapezoidal';	% 4 larger error than const
		%METHOD  = 'simpson';		% 0 not positive var
	end
	methods (Static)
		H = cdfS(h);
		cm = cmomentS(h,edge,order,v);
		[es eb] = entropyS(h,edge);
		ku = kurtosisS(h,edge,v);
		mu = meanS(h,edge,v);
		me = medianS(h,edge);
		mo = modeS(h,edge);
		m  = momentS(h,edge,order,v);
		q  = quantileS(h,edge,p);
		sk = skewnessS(h,edge,v);
		s  = stdS(h,edge,v);
		s2 = varS(h,edge,v);
		[wh e] = stairsS(h,edge);
		function centre = centreS(edge)
			centre = 0.5*(edge(1:end-1)+edge(2:end));
		end
		function [flag] = validS(h)
			s = sum(h,2);
			flag = isfinite(s) & (s > 0);
		end
	end
	methods
	function obj = Histogram(h,edge)
		obj.edge = rvec(edge);
		obj.h = h;
		obj.normalize();
	end
	function obj = normalize(obj)
            % normalise pdf
       	    obj.h(isnan(obj.h)) = 0;
            obj.h               = bsxfun(@times,obj.h,1./sum(obj.h,2));
	end
	function [flag] = valid(obj)
		s = sum(obj.h,2);
		flag = isfinite(s) & (s > 0);
	end

	% pseudo members
	function [centre obj] = centre(obj)
		centre = 0.5*(obj.edge(1:end-1)+obj.edge(2:end));
	end % centre
	end % methods
end % class Histogram

