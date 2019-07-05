% 2018-07-19 20:44:52.704290732 +0200
%
%% extension of the Theil-Senn regression to higher dimensions by
%% means of the Gauss-Seidel iteration
%
% TODO the algorithm converges slowly in certain situations, where the
%      paramters fluctuate around the optimum
classdef Theil_Multivariate < handle
	properties
		AbsTol  = 1e-7;
		RelTol  = 1e-3;	
		MaxIter = 1000;
		param   = [];
	end % properties
	methods
		function param = fit(obj,X,Y)
			c = zeros(size(X,2),1);

%note: these are n, not n^2 samples per iteration
%	-> this is basically gradient descent in one coordinate direction per iteration
%	-> this is more robust for small sample sizes (!)
%iterate
%for jth parameter
% for ith sample
%	-> determine dpj,i
%	0 = f0 + dfi_dxpj dpj,i
% end
% pj = pj + median(dpj,:)
%end

			% gauss seidel iteration
			% repeat until convergence
			iter = 0;
			while (1)
				iter = iter+1;
				C(:,iter) = c;
				cold = c;
				% y = sum si * xi
				res = Y-X*c;
				% keep all parameters fixed except one and optimise the p-th parameter
				for idx=1:size(X,2)
					c(idx) = obj.slope(X(:,idx), ...
							   res+X(:,idx)*c(idx));
				end % for idx
				dc = c-cold;
				if (    max(abs(dc)) < obj.AbsTol ...
				     || max(abs(dc./c)) < obj.RelTol ...
				   )
					printf('Iteration stopped after %d iterations, because iteration converged',iter);
					break;
				end % if convergence
				if ( iter > obj.MaxIter )
					printf('Iterations stopped, because maximum number of iterations was reached');
					figure(1)
					clf
					plot(C')
					pause

					break;
				end
			end % while 1
			if (iter > 10)
					figure(1)
					clf
					plot(C')
					pause
			end
			% constant
			c0    = hodges_lehmann_location(Y-X*c);
			param = [c0;c];
			obj.param = param;
			% TODO statistics
		end % fit
		function s = slope(obj,x,y)
			ss = bsxfun(@minus,y,y') ./ bsxfun(@minus,x,x');
			flag = triu(true(length(x)));
			s  = nanmedian(ss(flag));
		end
		function y = predict(obj,X)
			y = X*obj.param(2:end)+obj.param(1);
		end % predict
	end % methods
end % classdef Theil_Multivariate

