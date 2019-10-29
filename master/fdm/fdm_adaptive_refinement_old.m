% Sun Mar 25 01:26:19 MSK 2012
% Karl Kästner, Berlin

% requires similarity-symmetry transform to be undone before call
% X is cell array for each dimension
% X must contain boundary points l_x = l_v + 2
% n_d : size of hypercube, eg. [n_x] or [n_x n_y] or  [n_x, n_y, n_z], ...
% D2 : Laplacian difference matrix
% number k of new segments per element boundary
function [X_new err v4 mv4] = fdm_adaptive_refinement(X, n_d, v, D2, k) 
	% get fourth derivative by multiplying twice with the second order difference matrix
	v4 = D2*(D2*v);
	if (length(v4) > 2)
		v4(1) = 2*v4(2)-v4(3);
	end
	% interpolate derivatives to element centres
	% mv4 = max( [0; abs(v4(1:end))], [abs(v4(1:end)); 0] )
	% reshape into plane or cube
	if (length(n_d) > 1)
		v4 = reshape(v4, n_d);
	end
	% attach boundary zeros to v4
	v4 = padarray(v4,n_d.^0);		% TODO 4th derivative is not zero at the boundary
%	v4 = [2*v4(1)-v4(2); v4; 2*v4(end)-v4(end-1)];
	% interpolate fourth derivative into cell centres
	switch (length(n_d))
		case {1}  	
				X1 = X{1};
				mv4 = interp1(X1,v4, ...
					0.5*[X1(1:end-1)+X1(2:end)]);
		case {2}  
				X1 = X{1};
				X1 = X{2};
				mv4 = interp1(X1,X2,v4, ...
					0.5*[X1(1:end-1)+X1(2:end)], ...
					0.5*[X2(1:end-1)+X2(2:end)]);
		case {3}
			  	X1 = X{1};
				X2 = X{2};
				X3 = X{3};
				mv4 = interp1(X1,X2,X3,v4, ...
					0.5*[X1(1:end-1)+X1(2:end)], ...
					0.5*[X2(1:end-1)+X2(2:end)], ...
					0.5*[X3(1:end-1)+X3(2:end)]);
	end % switch
	% find the normed side length of each element
	h = zeros(size(mv4))/length(n_d);
	for ddx=1:length(n_d)
		n_d_ = n_d;
		n_d_(ddx) = 1;
		h = h + repmat(diff(X{ddx}), n_d_).^2;
	end % for ddx
	h = sqrt(h);
	% estimate the error contribution per element
	err = 1/12*h.^2.*mv4;
	err_abs = abs(err);
	for ddx=1:length(n_d)
		% find highest error in respective dimension
		err_max = max(err_abs,[],ddx);
		% mark cells for refinement
		% for k small elements replacing one old element, the error is expected to drop by 1/k
		% as it is proportional to the longest element side
		rdx = find(err_abs > 1/k^2*err_max);
		% generate k new points per marked element
		l = length(rdx);
		X_new_ = zeros((k-1)*l,1);
		X_ = X{ddx};
		for idx=1:k-1
			X_new_((idx-1)*l+1:(idx)*l,1) = 1/k*((k-idx)*X_(rdx) + idx*X_(rdx+1));
		end	
		% sort old and new points together
		X_new{ddx} = sort([X{ddx}; X_new_]);%.*(1 + 1e-4*randn(length(X_new_),1))]);
	end % for ddx
end % adaptive_refinement

