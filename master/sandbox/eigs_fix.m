
function varargout = eigs_fix(A,varargin)
	if (1 == size(A,1))
		if (nargout() > 1)
			varargout{1} = 1;
			varargout{2} = A/B;
			if (nargout() > 2)
			varagout{3} = 1; % flag for convergence
			end
		else			
			varargout{1} = A/B;
		end
	else
		varargout = eigs(A, varargin);
	end
end

