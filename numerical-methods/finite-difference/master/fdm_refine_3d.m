% Wed May 23 23:25:38 MSK 2012
% Karl Kästner, Berlin

function [X err_est err] = fdm_refine_3d(X, v, mode, cnorm)

	% refinement threshhold (compared to maximum error)
	p = 1/3;

	% factor to avoid "element" degeneration
	q = 1e-7;
	
	X1 = X{1};
	X2 = X{2};
	X3 = X{3};
	n1 = size(X1,1);
	n2 = size(X2,1);
	n3 = size(X3,1);

	if (0 == mode)
	% uniform refinement
		err_max = 1;
		err = [];
		X{1} = linspace(X1(1), X1(end), 2*n1)';
		X{2} = linspace(X2(1), X2(end), 2*n2)';
		X{3} = linspace(X3(1), X3(end), 2*n3)';
		return
	end

	I1 = speye(n1-2);
	I2 = speye(n2-2);
	I3 = speye(n3-2);

	% step widths
	H1 = diff(X1);
	H2 = diff(X2);
	H3 = diff(X3);

	% "element" volumes
	A = kron(kron(H1,H2),H3);

	% todo, correct bc in D4 and D3
	if (1 == mode)
		% estimate the discretisation error by higher order derivatives

		% compute 3rd and 4th order partial derivatives
		D3_1 = fdm_d_vargrid(X1, 3, 5);
		D3_1 = kron(kron(D3_1, I2), I3);
		v3_1 = D3_1*v;
		D4_1 = fdm_d_vargrid(X1, 4, 5);
		D4_1 = kron(kron(D4_1, I2), I3);
		v4_1 = D4_1*v;
	
		D3_2 = fdm_d_vargrid(X2, 3, 5);
		D3_2 = kron(kron(I1, D3_2), I3);
		v3_2 = D3_2*v;
		D4_2 = fdm_d_vargrid(X2, 4, 5);
		D4_2 = kron(kron(I1, D4_2), I3);
		v4_2 = D4_2*v;

		D3_3 = fdm_d_vargrid(X3, 3, 5);
		D3_3 = kron(kron(I1, I2), D3_3);
		v3_3 = D3_3*v;
		D4_3 = fdm_d_vargrid(X3, 4, 5);
		D4_3 = kron(kron(I1, I2), D4_3);
		v4_3 = D4_3*v;

		% step widths
		Hl1 = diff(X1(1:end-1));
		Hr1 = diff(X1(2:end));
		Hl2 = diff(X2(1:end-1));
		Hr2 = diff(X2(2:end));
		Hl3 = diff(X3(1:end-1));
		Hr3 = diff(X3(2:end));

		% estimate the error
		err1 = 1/3*abs(kron(kron(diag(sparse(Hl1-Hr1)), I2), I3)*v3_1) ...
			+ 1/12*abs(kron(kron(diag(sparse(Hl1.^2 - Hl1.*Hr1 + Hr1.^2)), I2), I3)*v4_1 );
		err2 = 1/3*abs(kron(kron(I1, diag(sparse(Hl2-Hr2))),I3)*v3_2) ...
			+ 1/12*abs(kron(kron(I1, diag(sparse(Hl2.^2 - Hl2.*Hr2 + Hr2.^2))),I3)*v4_2 );

		err3 = 1/3*abs(kron(kron(I1, I2), diag(sparse(Hl3-Hr3)))*v3_3) ...
			+ 1/12*abs(kron(kron(I1, I2), diag(sparse(Hl3.^2 - Hl3.*Hr3 + Hr3.^2)))*v4_3 );
	else
		% estimate the discretisation erro by
		% comparing three and five point kernels
		L1_3 = fdm_d_vargrid(X1,2,3);
		L1_5 = fdm_d_vargrid(X1,2,5);
		err1 = abs(kron(kron(L1_5 - L1_3, I2),I3)*v);
		L2_3 = fdm_d_vargrid(X2,2,3);
		L2_5 = fdm_d_vargrid(X2,2,5);
		err2 = abs(kron(kron(I1, L2_5 - L2_3),I3)*v);
		L3_3 = fdm_d_vargrid(X3,2,3);
		L3_5 = fdm_d_vargrid(X3,2,5);
		err3 = abs(kron(kron(I1, I2), L3_5 - L3_3)*v);
	end
	% interpolate error into "element" centres
	X1c = 0.5*(X1(1:end-1) + X1(2:end));
	X2c = 0.5*(X2(1:end-1) + X2(2:end));
	X3c = 0.5*(X3(1:end-1) + X3(2:end));
	[x y z] = meshgrid(X1(2:end-1), X2(2:end-1), X3(2:end-1));
	[xc yc zc] = meshgrid(X1c(2:end-1), X2c(2:end-1), X3c(2:end-1));
	err1 = reshape(err1,n3-2,n2-2,n1-2);
	err2 = reshape(err2,n3-2,n2-2,n1-2);
	err3 = reshape(err3,n3-2,n2-2,n1-2);
	err1 = shiftdim(err1,1);
	err2 = shiftdim(err2,1);
	err3 = shiftdim(err3,1);
%clf; subplot(1,2,1)
%slice(x, y, z, err1, 1, 1, 1);
	err1 = interp3(x,y,z, err1, xc,yc,zc);
%subplot(1,2,2)
%slice(xc, yc, zc, err1, 1, 1, 1);
%pause
	err2 = interp3(x,y,z, err2, xc,yc,zc);
	err3 = interp3(x,y,z, err3, xc,yc,zc);
 	% extrapolate for boundary elements
	err1 = padarray(err1,[1 1 1],'replicate');
	err2 = padarray(err2,[1 1 1],'replicate');
	err3 = padarray(err3,[1 1 1],'replicate');
	% weight error by "element" area (L1-norm)
	A = reshape(A,n3-1,n2-1,n1-1);
	A = shiftdim(A,1);
	switch (cnorm)
		case{1} % L1-norm
			% weight by area
			err1 = A.*err1;
			err2 = A.*err2;
			err3 = A.*err3;
			% combined error of all coordinate axis (for plot)
			err = err1 + err2 + err3;
			% summed error per dimension
			err1 = squeeze(sum(sum(err1,3),1))';
			err2 = squeeze(sum(sum(err2,3),2));
			err3 = squeeze(sum(sum(err3,2),1));
		otherwise % L_inf norm
			% combined error of all coordinate axis (for plot)
			err = err1 + err2 + err3;

			% maximum error per row / column
			err1 = squeeze(max(max(err1,[],3),[],1))';
			err2 = squeeze(max(max(err2,[],3),[],2));
			err3 = squeeze(max(max(err3,[],2),[],1));
	end

	% estimated error ( convergence criteria)
		%area = sum(sum(sum(A)));
	err_est = sum(sum(sum(err)));
	
	% threshold for refinement
	err_max = max([max(err1) max(err2) max(err3)]);
	%err3 = max(squeeze(max(err3,[],2))',[],1) % same
%clf; %subplot(1,2,1)
%[xc yc zc] = meshgrid(X1c, X2c, X3c);
%subplot(2,2,1)
%slice(xc, yc, zc, err1, 0, 0, 0);
%subplot(2,2,2)
%slice(xc, yc, zc, err2, 0, 0, 0);
%subplot(2,2,3)
%slice(xc, yc, zc, err3, 0, 0, 0);
%subplot(2,2,4)
%clf
%size(err)
%slice(xc, yc, zc, log10(err), 0, 0, 0);
%pause

	% minimum grid spaces
	h1min = min(H1);
	h2min = min(H2);
	h3min = min(H3);

	% refine in x-direction
	X_ = zeros(n1,1);
	n_ = 0;
	for idx=1:n1-1
		if (err1(idx,1) > p*err_max || q*H1(idx) > min(h2min,h3min));
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
		if (err2(idx,1) > p*err_max || q*H2(idx) > min(h1min,h3min));
			X_(n_+1,1) = 1/3*X2(idx,1) + 2/3*X2(idx+1,1);
			X_(n_+2,1) = 2/3*X2(idx,1) + 1/3*X2(idx+1,1);
			n_ = n_+2;
		end % if
	end % for idx
	X2 = sort([X2; X_(1:n_)]);

	% refine in z-direction
	X_ = zeros(n3,1);
	n_ = 0;
	for idx=1:n3-1
		if (err3(idx,1) > p*err_max || q*H3(idx) > min(h1min,h2min));
			X_(n_+1,1) = 1/3*X3(idx,1) + 2/3*X3(idx+1,1);
			X_(n_+2,1) = 2/3*X3(idx,1) + 1/3*X3(idx+1,1);
			n_ = n_+2;
		end % if
	end % for idx
	X3 = sort([X3; X_(1:n_)]);

	% stack back
	X{1} = X1;
	X{2} = X2;
	X{3} = X3;
end % fdm_refine_3d

