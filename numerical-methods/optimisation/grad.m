% Mon Aug 25 18:13:15 CEST 2014
% Karl Kastner, Berlin
%
%% numerical gradient
%
function [g f0] = grad(func,x,h,mode,parallel)
    if (nargin()<3 || isempty(h))
		h = abs(x).'*sqrt(eps); % class(x);
		h = max(h,sqrt(eps));
		%h = max(h,1e3*sqrt(eps));
    else
        h = rvec(h);
	end
	if (1 == length(h))
		h = h*ones(1,length(x));
	end

	if (nargin()<4 || isempty(mode))
		mode = '';
	end
	if (nargin() < 5)
		parallel = false;
	end
	switch (mode)
	case {'one-sided'}
		f0 = feval(func,x);
		fr = zeros(length(f0),length(x),class(f0));
		if (parallel)
			parfor idx=1:length(x)
				x_        = x;
				x_(idx)   = x_(idx)+h(idx);
				fr(idx,1) = feval(func,x_);
			end
		else
			for idx=1:length(x)
				x_        = x;
				x_(idx)   = x_(idx)+h(idx);
				fr(:,idx) = feval(func,x_);
			end
		end
		g = bsxfun(@times,(1./h),bsxfun(@minus,fr,f0));
	otherwise % centred
		% TODO allocate memory
		fl = []; %zeros(length(f),length(x)); %,class(f0));
		fr = [];% zeros(length(f),length(x)); %,class(f0));
		for idx=1:length(x)
			xr = x;
			xr(idx)   = xr(idx)+h(idx);
			xl = x;
			xl(idx)   = xl(idx)-h(idx);
			fr(:,idx) = feval(func,xr);
			fl(:,idx) = feval(func,xl);
		end
		f0 = 0.5*(fl+fr);
		g = bsxfun(@times,(0.5./h),(fr-fl));
	end
end % function grad

