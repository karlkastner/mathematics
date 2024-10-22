% Thu  3 Oct 10:08:54 CEST 2024
% Karl Kastner, Berlin
%% check convergence requirements of coarsened matrix
function check(obj)
	n = obj.n;
	for k=1:length(obj.s)
		J = obj.s(k).diagonals;
		D = abs(J(:,5));
		R = sum(abs(J(:,1:4)),2);

		% check diagonal dominance condition for convergence
		% note that the iteration can also converge when this condition
		% is not met as long J is SPD
		maxratio = max(R./D);
		% check that column sum is zero
		% note that this only holds at interior points
		% it holds only for boundary points with Neumann or circular boundary conditions
		% TODO, these are rows, not columns
		J = diags2mat(J,n);
		abssumcol  = abs(sum(J,1)-1);
		sumabscol = sum(abs(J),1)-1;
		maxrelcolsum = max(abssumcol./sumabscol);
		eigs(@(x) bicgstabl(J,x,100),prod(n),2,'largestabs')
		eigs(@(x) bicgstabl(J,x,100),prod(n),2,'largestreal')
		printf('level %g n %4g %4g max |R|/|D| = %1.1e, max |sum J-I| %e max |sum J-I|/sum|J-I| %e\n',k,n,maxratio,max(abssumcol),maxrelcolsum);
		n = n/2;
	end
end


