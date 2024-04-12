-> so dx and dx^2 are correctly resampled
L=15; n = L+1; A = derivative_matrix_1_1d(n,L,2,'circular'); A = full(A), D=downsampling_matrix_1d(n,'mg'); U=upsampling_matrix_1d(n); full(A), D*A*U

