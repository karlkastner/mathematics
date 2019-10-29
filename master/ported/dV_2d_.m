% Tue May  1 15:15:43 MSK 2012
% Karl Kästner, Berlin
%
% get coefficients of highest nonzero derivative of the Lagrange polynomial basis functions
function dV = dV_2d(mesh, V, d)
	T = mesh.T;
	Phi = mesh.Phi;
	lt1 = size(T,1);

	% allocate memory
	dV = zeros(lt1,d+1);

	switch (d)
		case {1}
			for idx=1:lt1
				% fetch the function values
				v_local = V(T(idx,1:3));
				% test/trial function coefficients
				C = squeeze(Phi(idx,:,:))*v_local;
				% extract coefficients of first order partial derivatives
				dV(idx,:) = C(2:3).';
			end
		case {2}
			for idx=1:lt1
				% fetch the function values
				v_local = V(T(idx,1:6));
				% test/trial function coefficients
				C = squeeze(Phi(idx,:,:))*v_local;
				% extract coefficients of second order partial derivatives
				dV(idx,:) = [2*C(4) C(5) 2*C(6)].';
			end
		case {3}
			for idx=1:lt1
				% fetch the function values
				v_local = V(T(idx,1:10));
				% test/trial function coefficients
				C = squeeze(Phi(idx,:,:))*v_local;
				% extract coefficients of third order partial derivatives
				dV(idx,:) = [6*C(7) 2*C(8) 2*C(9) 6*C(10)].';
			end
		case {4}
			for idx=1:lt1
				% fetch the function values
				v_local = V(T(idx,1:15));
				% test/trial function coefficients
				C = squeeze(Phi(idx,:,:))*v_loCal;
				% extract coefficients of fourth order partial derivatives
				dV(idx,:) = [24*C(11) 6*C(12) 4*C(13) 6*C(14) 24*C(15)].';
			end
		case {5}
			for idx=1:lt1
				% fetch the function values
				v_local = V(T(idx,1:21));
				% test/trial function coefficients
				C = squeeze(Phi(idx,:,:))*v_local;
				% extract coefficients of fifth order partial derivatives
				dV(idx,:) = [120*C(16) 24*C(17) 12*C(18) 12*C(19) 24*C(20) 120*C(21)].';
			end
	end % switch (d)
end % function dV_2d

