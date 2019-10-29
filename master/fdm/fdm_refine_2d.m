% Mon May 21 23:26:52 MSK 2012
% Karl Kästner, Berlin

function [X err_est err] = fdm_refine_2d(X, v, mode, cnorm)

	% refinement threshhold (compared to maximum error)
	p = 1/3;
	if (cnorm < inf)
		p = p^cnorm;
	end

	% factor to avoid degeneration
	q = 1e-7;

	X1 = X{1};
	X2 = X{2};
	n1 = size(X1,1);
	n2 = size(X2,1);

	if (0 == mode)
	% uniform refinement
		err_est = 1;
		err = [];
		X{1} = linspace(X1(1), X1(end), 2*n1)';
		X{2} = linspace(X2(1), X2(end), 2*n2)';
		return
	end

	% identity matrices
	I2 = speye(n2-2);
	I1 = speye(n1-2);
	
	% step widths
	H1 = diff(X1);
	H2 = diff(X2);
	% "element" areas
	A = kron(H1,H2);

	% todo, correct bc in D4 and D3
	if (1 == mode)
		% estimate the discretisation error by higher order derivatives

		% compute 3rd and 4th order partial derivatives
		D3_1 = fdm_d_vargrid(X1, 3, 5);
		D3_1 = kron(D3_1, I2);
		v3_1 = D3_1*v;
		D4_1 = fdm_d_vargrid(X1, 4, 5);
		D4_1 = kron(D4_1, I2);
		v4_1 = D4_1*v;
	
		D3_2 = fdm_d_vargrid(X2, 3, 5);
		D3_2 = kron(I1, D3_2);
		v3_2 = D3_2*v;
		D4_2 = fdm_d_vargrid(X2, 4, 5);
		D4_2 = kron(I1, D4_2);
		v4_2 = D4_2*v;

		% step widths
		Hl1 = diff(X1(1:end-1));
		Hr1 = diff(X1(2:end));
		Hl2 = diff(X2(1:end-1));
		Hr2 = diff(X2(2:end));

		% estimate the error
		err1 = 1/3*abs(kron(diag(sparse(Hl1-Hr1)), I2)*v3_1) ...
			+ 1/12*abs(kron(diag(sparse(Hl1.^2 - Hl1.*Hr1 + Hr1.^2)), I2)*v4_1 );
		err2 = 1/3*abs(kron(I1, diag(sparse(Hl2-Hr2)))*v3_2) ...
			+ 1/12*abs(kron(I1, diag(sparse(Hl2.^2 - Hl2.*Hr2 + Hr2.^2)))*v4_2 );
	else
		% estimate the discretisation erro by
		% comparing three and five point kernels
		L1_3 = fdm_d_vargrid(X1,2,3);
		L1_5 = fdm_d_vargrid(X1,2,5);
		err1 = abs(kron(L1_5 - L1_3, I2)*v);
		L2_3 = fdm_d_vargrid(X2,2,3);
		L2_5 = fdm_d_vargrid(X2,2,5);
		err2 = abs(kron(I1, L2_5 - L2_3)*v);
	end

	% interpolate error into "element" centres
	X1c = 0.5*(X1(1:end-1) + X1(2:end));
	X2c = 0.5*(X2(1:end-1) + X2(2:end));
	err1 = reshape(err1,n2-2,n1-2);
	err1 = interp2(X1(2:end-1), X2(2:end-1)', err1, X1c(2:end-1), X2c(2:end-1)');
	err2 = reshape(err2,n2-2,n1-2);
	err2 = interp2(X1(2:end-1), X2(2:end-1)', err2, X1c(2:end-1), X2c(2:end-1)');
 	% extrapolate for boundary elements
	err1 = padarray(err1,[1 1],'replicate');
	err2 = padarray(err2,[1 1],'replicate');

	% estimated error
	A = reshape(A,n2-1,n1-1);
	area = sum(sum(A));
	% A = (1/area)*A;
	switch (cnorm)
		case {1}
			% weight error by "element" area (L1-norm)
			% (error contribution to global error)
			err1 = A.*abs(err1);
			err2 = A.*abs(err2);

			% combined error of all coordinate axis (for plot)
			err = err1 + err2;

			% convergence criteria
			err_est = sum(sum(err));

			% summed error for each dimension
			err1 = sum(err1,1)';
			err2 = sum(err2,2);
		case {2}
			% weight squared error by "element" area (L2-norm)
			% (error contribution to global error)
			err1 = A.*err1.^2;
			err2 = A.*err2.^2;

			% combined error of all coordinate axis (for plot)
			err = err1 + err2;

			% convergence criteria
			err_est = sqrt(sum(sum(err))*area);

			% summed error for each dimension
			err1 = sum(err1,1)';
			err2 = sum(err2,2);
		otherwise % choose infinity norm
			% do not weight the error
			err1 = abs(err1);
			err2 = abs(err2);

			% combined error of all coordinate axis (for plot)
			err = err1 + err2;

			% convergence criteria
			err_est = max(max(err));

			% maximum error per row / column
			err1 = max(err1,[],1)';
			err2 = max(err2,[],2);
	end % switch cnorm

	% threshold for refinement is p*err_max
	err_max = max(max(err1), max(err2));

	% minimum step widths
	h1min = min(H1);
	h2min = min(H2);

	% refine in x-direction
	X_ = zeros(n1,1);
	n_ = 0;
	for idx=1:n1-1
		if (err1(idx,1) > p*err_max || q*H1(idx) > h2min);
			X_(n_+1,1) = 1/3*X1(idx,1) + 2/3*X1(idx+1,1);
			X_(n_+2,1) = 2/3*X1(idx,1) + 1/3*X1(idx+1,1);
			n_ = n_+2;
		end % if
	end % for idx
	X1 = sort([X1; X_(1:n_)]);


	% refine in y-direction
	X_ = zeros(n2,1);
	n_ = 0;
	for idx=1:n2-1
		h=1;
		if (err2(idx,1) > p*err_max || q*H2(idx) > h1min);
			X_(n_+1,1) = 1/3*X2(idx,1) + 2/3*X2(idx+1,1);
			X_(n_+2,1) = 2/3*X2(idx,1) + 1/3*X2(idx+1,1);
			n_ = n_+2;
		end % if
	end % for idx
	X2 = sort([X2; X_(1:n_)]);

	% smoothing
%	n = length(X1)+1;
%	M=full(spdiags( [ [0.5*ones(n-3,1); 0; 0] [1; zeros(n-3,1); 1] [0; 0; 0.5*ones(n-3,1)]], -1:1, n-1,n-1));
%	X1 = M*X1;
%	X1 = M*X1;
%	X1 = M*X1;
%	n = length(X2)+1;
%	M=full(spdiags( [ [0.5*ones(n-3,1); 0; 0] [1; zeros(n-3,1); 1] [0; 0; 0.5*ones(n-3,1)]], -1:1, n-1,n-1));
%	X2 = M*X2;
%	X2 = M*X2;
%	X2 = M*X2;

	% stack back
	X{1} = X1;
	X{2} = X2;

end % fdm_refine_2d

