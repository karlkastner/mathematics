% triangle scheme - different schemes for different order of accuracy
% returns phi and d-phi value at coordinate for triangle k in T
[phi dphi] = phi(k, P, T, x, y)
% integration scheme - different schemes depending on order of accuracy
% returns partial cell of value A_ij, pi,pj \in 1..3
% f   : coefficient function
% phi : trial function
[v] = int(nt, P, T, pi, pj, phi, f)
%assemble scheme
%returns matrix
	for each triangle k
		for each of the three test functions	i
			for each of the three test functions	j
				A_ij = int(k, P, T, i, j, phi, f)

