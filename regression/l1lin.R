# 2014-08-25 22:43:03 +0200
# Karl Kastner, Berlin
#
## minimise ||Ax - b||_L1
#
l1lin <- function (A,b)
{
	# number of parameters
	nx = dim(A)[1];
	# number of samples
	nb = dim(A)[2];

	# matrix casting the nx parameter minimisation problem
	# into an nx+nb parameter linear programming problem
	I  = diag(1,nx);
	Z  = matrix(0,nx,nb);
	AA = rbind(rbind(cbind(-A, -I),
	                 cbind( A, -I)),
	                 cbind( Z,  I));
	bb = rbind(rbind(-b,b),matrix(0,nb,1));

	# weights for optimisation	
	f = rbind(matrix(1,np,1),
	          matrix(0,nx,1));
	# prepare and pass arguments to linear programming solver
	lprec <- make.lp(nrow=dim(AA)[1], ncol=np+nx)
	set.objfn(lprec, f)
	for (i in 1:dim(AA)[1])
	{
		add.constraint(lprec, AA[i,], "<=", bb[i]);
	}
	# minimise
	solve(lprec)
	return (lprec);
} # function l1lin

