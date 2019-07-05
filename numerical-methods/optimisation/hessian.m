% Fr 5. Feb 15:41:30 CET 2016
% Karl Kastner, Berlin
%
%% numerical hessian
% function [H g f0 X F] = hessian(fun,x0,h)
%
function [H, g, f0, X, F] = hessian(fun,x0,h)
	parflag = false; %true;
	n = length(x0);
	if (nargin() < 3 || isempty(h))
		% note: h < eps^0.5 for consistent hessian !!!
		% h for gradient can be eps^0.5
		h = max(abs(x0)*eps.^0.25,eps.^0.25);
	end
	if (isscalar(h))
		h = h*ones(n,1);
	end
	Eh  = diag(sparse(h));

	fl = zeros(n,1); 
	fr = zeros(n,1); 
	% TODO only two matrices are necessary, as only the upper triangle is used
	fll = zeros(n,n);
	flr = zeros(n,n);
	frl = zeros(n,n);
	frr = zeros(n,n);
if (~parflag)
	% function at origin
	f0 = fun(x0);

	for idx=1:n
		% left function value
		fl(idx) = fun(x0-Eh(:,idx));
		% right function value
		fr(idx) = fun(x0+Eh(:,idx));		
		% off diagonal values
		for jdx=idx+1:n
			fll(idx,jdx) = fun(x0-Eh(:,idx)-Eh(:,jdx));
			flr(idx,jdx) = fun(x0-Eh(:,idx)+Eh(:,jdx));
			frl(idx,jdx) = fun(x0+Eh(:,idx)-Eh(:,jdx));
			frr(idx,jdx) = fun(x0+Eh(:,idx)+Eh(:,jdx));	
		end
	end
else
	parfor pdx=1:n*n
		% function value at origin                                                    
		if (pdx <= 1)
			f0(pdx) = fun(x0); 
		end
		% for gradient and diagonal of the hessian
		if (pdx <= n)
			% left function value
			fl(pdx) = fun(x0-Eh(:,pdx));
			% right function value
			fr(pdx) = fun(x0+Eh(:,pdx));		
		end % if
		% for off-diagonal of the hessian
		[idx jdx] = ind2sub([n n],pdx);
		if (jdx > idx)
			fll(pdx) = fun(x0-Eh(:,idx)-Eh(:,jdx));
			flr(pdx) = fun(x0-Eh(:,idx)+Eh(:,jdx));
			frl(pdx) = fun(x0+Eh(:,idx)-Eh(:,jdx));
			frr(pdx) = fun(x0+Eh(:,idx)+Eh(:,jdx));
		end % if
		%end % else of if
	end % for pdx
end

	% gradient
	g = 0.5*(fr-fl)./h;

	% diagonal and upper triangle of the Hessian
	H =     diag((fl-2*f0+fr)./(h.*h));
	H = H + 0.25*((frr-frl) - (flr-fll))./(h*h');
	% lower triangle of the hessian
	H = H + triu(H,+1)';

	if (nargout() > 3)
		fdx = triu(true(n),+1);
		% concatenate x-values
		Xl = [];
		Xr = [];
		Xll = [];
		Xlr = [];
		Xrl = [];
		Xrr = [];
		for idx=1:n
			Xl(:,idx) = x0-Eh(:,idx);
			Xr(:,idx) = x0+Eh(:,idx);
			for jdx=idx+1:n
				Xll = [Xll,x0-Eh(:,idx)-Eh(:,jdx)];
				Xlr = [Xlr,x0-Eh(:,idx)+Eh(:,jdx)];
				Xrl = [Xrl,x0+Eh(:,idx)-Eh(:,jdx)];
				Xrr = [Xrr,x0+Eh(:,idx)+Eh(:,jdx)];
			end
		end
		X = [x0,Xl,Xr,Xll,Xlr,Xrl,Xrr];
		% concatenate function values
		F = [f0;fl;fr;fll(fdx);flr(fdx);frl(fdx);frr(fdx)];
	end
end % hess

