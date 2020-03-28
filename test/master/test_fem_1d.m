% Wed Jul 11 21:32:37 MSK 2012
% Karl Kästner, Berlin

function test_fem_1d()

N = 2.^(0:6);
func = [];

e_true = -pi^2;

for idx=1:length(N)
	n = N(idx);
	[P T BC] = mesh_1d_uniform(n+2,1);

	int = @int_1d_trapezoidal;
	A = assemble_1d_dphi_dphi(P, T, int, func);
	B = assemble_1d_phi_phi(P, T, int, func);
	[A B] = boundary_1d(A,B,BC,0);	
	E(idx,1) = eigs(A,B,1,'SM');
	L{1} = 'linear, trapezoidal';

	int = @int_1d_gauss_1;
	A = assemble_1d_dphi_dphi(P, T, int, func);
	B = assemble_1d_phi_phi(P, T, int, func);
	[A B] = boundary_1d(A,B,BC,0);	
	E(idx,2) = eigs(A,B,1,'SM');
	L{2} = 'linear, gauss';

	int = @int_1d_gauss_2;
	[P_ T_ BC_] = promote_1d_2_3(P, T, BC);
	A = assemble_1d_dphi_dphi(P_, T_, int, func);
	B = assemble_1d_phi_phi(P_, T_, int, func);
	[A B] = boundary_1d(A,B,BC_,1);	
	E(idx,3) = eigs(A,B,1,'SM');
	L{3} = 'quadratic, gauss';
	
	int = @int_1d_gauss_3;
	[P_ T_ BC_] = promote_1d_2_4(P, T, BC);
	A = assemble_1d_dphi_dphi(P_, T_, int, func);
	B = assemble_1d_phi_phi(P_, T_, int, func);
	[A B] = boundary_1d(A,B,BC_,1);
	E(idx,4) = eigs(A,B,1,'SM');
	L{4} = 'cubic, gauss';
	
	int = @int_1d_gauss_4;
	[P_ T_ BC_] = promote_1d_2_5(P, T, BC);
	A = assemble_1d_dphi_dphi(P_, T_, int, func);
	B = assemble_1d_phi_phi(P_, T_, int, func);
	[A B] = boundary_1d(A,B,BC_,1);
	E(idx,5) = eigs(A,B,1,'SM');
	L{5} = 'quartic, gauss';
	
	int = @int_1d_gauss_5;
	[P_ T_ BC_] = promote_1d_2_6(P, T, BC);
	A = assemble_1d_dphi_dphi(P_, T_, int, func);
	B = assemble_1d_phi_phi(P_, T_, int, func);
	[A B] = boundary_1d(A,B,BC_,1);
	E(idx,6) = eigs(A,B,1,'SM');
	L{6} = 'quintic, gauss';
end % for idx
E
E./e_true-1
Err = E - e_true;

loglog(N,abs(Err));
legend(L);
ylabel('|\lambda_1^* - \lambda_1^h|')
xlabel('number of grid points n')

end % test_fem_1d

