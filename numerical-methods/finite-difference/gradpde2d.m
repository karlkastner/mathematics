% Thu 17 May 16:51:21 CEST 2018
%
%% objective function gradiend on two dimensional regular grid
%% numeric gradient for non-linear least squares optimisation
%% of a PDE on a rectangular grid
%% x_* = min(f(x))
%% f = (v(x) - v(x_*))^2 = f(x) + A dx + O(dx^2)
%% a_ij = df_i/dx_j
%% 
%
function [f0,A] = gradpde2d(fun,x0,n)
	% function value at initial point
	f0 = fun(x0);

	% step for gradient computation
	% this is a hessian matrix, so sqrt(eps) is too small
	dx  = (eps^0.125)*max(1.0,abs(x0)).^0;
	%dx_ = zeros(n(1)+2,n(2)+2);
	%dx_(2:n(1)+1,2:n(2)+1) = dx;

	% buffered sparse matrix
	nn = prod(n+2);
	A  = sparse([],[],[],nn,nn,5*nn);
	nn = [nn,nn];

	df = zeros(n(1)+2,n(2)+2);
	df_ = 0;

	% a perturbation of a value influences approximately only the immediate
	% eight neighbours, so every third function value can be perturbed
	% to evaluate the perturbed function value in parallel
	% for the leading three columns
	s = 5;
	for col=1:s
	 % for the leading three rows
	 for row=1:s
		% group index
		r  = row:s:n(1);
		c  = (col:s:n(2))';
		id  = flat(bsxfun(@plus,r,      (c-1)*(n(1))));
		% buffered index
		id_ = flat(bsxfun(@plus,r+1,(c)*(n(1)+2)));
		% perturb x
		xp     = x0;
		xp(id) = xp(id) + dx(id);
		% function value
		fp = fun(xp)
		% difference in function value
		df = fp-f0;
		df_(id) = df(id);
		df(2:n(1)+1,2:n(2)+1) = reshape(df,n);
		
		%df(:,1) = 1;
		% entries of gradient matrix
		% self
		A(sub2ind(nn,id_,id_)) = df(id_)./dx(id);
		% left
		jd  = id -n(1);
		jd_ = id_-n(1)-2;
		A(sub2ind(nn,jd_,id_)) = df(jd_)./dx(id);
		% right
		jd  = id +n(1);
		jd_ = id_+n(1)+2;
		A(sub2ind(nn,jd_,id_)) = df(jd_)./dx(id);
		% up
		jd  = id -1;
		jd_ = id_-1;
		A(sub2ind(nn,jd_,id_)) = df(jd_)./dx(id);
		% down
		jd  = id +1;
		jd_ = id_+1;
		A(sub2ind(nn,jd_,id_)) = df(jd_)./dx(id);
         end % row
        end % col

%	figure();
%	clf()
%	imagesc(reshape(df_,n))
%	'honk'
%pause

	% remove buffer
	inner = false(n+2);
	inner(2:end-1,2:end-1) = true;
	inner = flat(inner);
	A = A(inner,inner);
%	A = A(2:nn-1,2:nn-1);

	% restore symmetry
	% A = 1/2*(A+A');
end % gradpde2d

