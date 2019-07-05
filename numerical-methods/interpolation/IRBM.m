% Sun Jul 20 16:09:13 WIB 2014
% Karl Kastner, Berlin
%
%% interpolate by the radial basis function method
% interpolate value from radial basis functions
% if V is multidimensional, no unique minimisation radius exists!!! -> hyperellipsoid
function [Vi Err R2] = IRBM(X,V,Xi,n_order,r2max)
	addpath([ROOTFOLDER,fielsep(),'../master/thesis/src/sandbox/ported/']);
	
	% allocate memory
	Vi  = NaN(size(Xi,1),size(V,2),class(V));
	Err = NaN(size(Xi,1),size(V,2),class(V));

	n_iter = 20;
	% TODO this implementation is for 2D only so far
	n_dim = size(Xi,2);

	% initial radius
%	xm = mean(Xi,1);
	if (nargin() < 5)
		xs = var(Xi,[],1);
		r2max = 4*geomean(xs)/size(Xi,1).^2;
	end
%	r2max  = max( (Xi(:,1) - xm(1)).^2 + (Xi(:,2) - xm(2)).^2 );
	R2  = repmat(r2max,size(Xi,1),1);
	
	% TODO convolve with cosine filter
	% for each target point
	t0 = tic();
	for idx=1:size(Xi,1)
		if (toc(t0) > 10)
			fprintf(1,'Progress IRBM: %d%%\n',round(100*idx/size(Xi,1)));
			t0 = tic();
		end
		% distance of source points to target points
		D1  = X(:,1) - Xi(idx,1);
		D2  = X(:,2) - Xi(idx,2);
		R2i = D1.^2 + D2.^2;
		for jdx=1:n_iter
			% find all points within the radius
			fdx = find(R2i < R2(idx,jdx));
			% TODO, minimum number of points
			if (length(fdx) < 6)
				Err(idx,jdx)  = NaN;
				R2(idx,jdx+1) = 2*R2(idx,jdx);
				if (jdx > 1)
					Vi(idx,jdx) = Vi(idx,jdx-1);
				else
					Vi(idx,jdx) = NaN;
				end
			else	
			% cosine filter
			% TODO this does not account for correlation
			% to assure inifnite derivatives and compact support
			%r = sqrt(R2(idx,jdx));
			%W = cos(D1(fdx)/r*pi).*cos(D2(fdx)/r*pi);
			w = 1 - sin(0.5*pi*sqrt(R2i(fdx)./R2(idx,jdx))).^2;
			w = w/sum(w);
			W = diag(sparse(double(w)));
			% cos-weight makes it unstable why?
			W = speye(size(W));
			% vandermonde matrices
			A  = vander_2d([D1(fdx) D2(fdx)],n_order);
			A_ = vander_2d([D1(fdx) D2(fdx)],n_order+1);
			if (rcond(A_'*W*A_) < 1e-6)
				Err(idx,jdx)  = NaN;
				R2(idx,jdx+1) = 2*R2(idx,jdx);
				if (jdx > 1)
					Vi(idx,jdx) = Vi(idx,jdx-1);
				else
					Vi(idx,jdx) = NaN;
				end
				continue
			end
			% calculate the interpolation coefficients
			p = (A'*W*A) \ (A'*W*V(fdx));
			p_ = (A_'*W*A_) \ (A_'*W*V(fdx));
%			vi = func(X(fdx,:),V(fdx,:),Xi(idx,2),W);\
			% evaluate polynomial at targe point (origin)
			Vi(idx,jdx) = p(1);
			Vi_(idx,jdx) = p_(1);
			% evaluate polynomial at source points
			vi = A*p;
			vi_ = A_*p_;
			% interpolation error
		%	err_i = Vi__ - Vi_;
		%	err_i = abs(p(1) - p_(1));
			% TODO, weights have to play a role here
			err_i = std(vi - vi_);
			% error due to noise
			err_n = std(vi_ - V(fdx))/sqrt(length(fdx));
		%	err_n = Vi_ - V(fdx,:);
			% the error is the sum of the interpolation error and noise
			Err_i(idx,jdx) = err_i; %'*err_i;
			Err_n(idx,jdx) = err_n; %'*err_n;
			% under assumption of no correlation
			%Err(idx,jdx) = sqrt(err_i'*W*err_i + err_n'*W*err_n);
			Err(idx,jdx) = sqrt(err_i'*err_i + err_n'*err_n);
			% optimal number of support points
			% noise contribution to the error becomes less for increasing radius,
			% whilst interpolation contribution becomes more for increasing radius
			p = 0.5;
			R2(idx,jdx+1) = R2(idx,jdx)*min(2,max(0.5,p*err_n/err_i + (1-p)));
		%	n_opt = err_i'*err_i/(err_n'*err_n);
			%n_opt = err_n'*err_n/(err_i'*err_i);
			% compute optimal radius
			% TODO, this is nD but a linear increase of points is assumed
		%	n = length(fdx);
		%	p = 0.5;
		%	r2_opt = p*R2(idx,jdx)*n*n/(n_opt*n_opt) + (1-p)*R2(idx,jdx);
			% limit step
		%	R2(idx,jdx+1) = max(0.5*R2(idx,jdx),min(2*R2(idx,jdx),r2_opt));
			end % enough points
		end % for jdx
	end % for idx
%	'Honk'
%	subplot(4,1,1)
%	plot(sqrt(R2))
%	subplot(4,1,2)
%	plot(Err_i)
%	subplot(4,1,3)
%	plot(Err_n)
%	subplot(4,1,4)
%	plot(Vi)
%	ylim([0 20])
%	Vi
%	pause
end % IRBM

