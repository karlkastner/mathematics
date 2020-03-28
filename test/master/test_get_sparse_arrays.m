	clear
	n = 4;

	B = sparse(rand(n));
	B=B+B';

	[IR IC BD B_dat] = generate(B,[],0),

	[N_ IC_ IR_ B_dat_] = get_sparse_arrays(tril(B,-1));
	N_     = N_'
	IR_    = IR_'+1
	IC_    = diff(IC_)'
	AD_ = full(diag(B))'
	B_dat_ = B_dat

%	N_ = free_matrix(N_);
%	IR_ = free_matrix(IR_);
%	IC_ = free_matrix(IC_);
%	B_mat_ = free_matrix(B_dat_);

