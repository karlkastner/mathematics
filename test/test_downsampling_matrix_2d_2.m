% Tue 17 Sep 08:40:59 CEST 2024
% Karl Kastner, Berlin

% Results:
% for the fw-prolongation, the advection part converges to the unstable central scheme
% for the hw-prolongation, the advection part vanishes
% -> scheme required that pays attention to advection and diffusion


 rng(0);

 %P = downsampling_matrix_2d_hw(n)/8;
mode_C = {'fw','hw'}

for idx=1:length(mode_C)
 n=[8,8];
mode = mode_C{idx}

 [R1_,R2t_,R1,R2t] = downsampling_matrix_2d(n,mode);
 [R1__,R2t__] = downsampling_matrix_2d(n,'fw');


disp('Testing downsampling of data')
x = rand(n);
switch (mode)
case {'fw'}
	x_d  = reshape(R1_*x(:),n(1)/2,n(2)/2);
	x_du = reshape(R2t_*x_d(:),n(1),n(2));
	x_ = (left(x,1) + 2*x + right(x,1))/4;
	x_ = x_(:,1:2:end);
	x_ = (up(x_,1) + 2*x_ + down(x_,1))/4;
	x_d_ = x_(1:2:end,:);
	% same as: x_ = downsample_2d(x)
	x_du_ = upsample_2d(x_d_);
	
case {'hw'}
	x_d  = reshape((1/2*R1_+R2t_'/4)*x(:),n(1)/2,n(2)/2);
	x_du = reshape(4*(1/2*R1_+R2t_'/4)'*x_d(:),n(1),n(2));

	x_ = (4*x+up(x,1)+down(x,1)+left(x,1)+right(x,1))/8;
	x_ = x_(1:2:end,1:2:end);
	x_d_ = x_;
	x_du_ = upsample_2d(x_d_);
end
rms(x_d-x_d_,'all')
rms(x_du-x_du_,'all')
x_du
x_du_
x_du - x_du_
disp('Testing downsampling of operators')
disp('asymptotic kernel test');

	m =10;
	n=2^m*[1,1];
	[Dx,Dy,D2x,Dxy,D2y] = derivative_matrix_2d(n,n-1,2,{'circular','circular'});
	[Dxn] = derivative_matrix_2d(n,n-1,-1,{'circular','circular'});
	
	D2=0.*(D2x+D2y)+Dxn;
	K  = extract_kernel(D2,n,[2,2])
	for idx=1:(m-2);
		[R1_,R2t_,R1,R2t] = downsampling_matrix_2d(n,mode);
		D2= (R1_*D2*R2t_);
		n = n/2;
	end;
	K_=extract_kernel(D2,n,[2,2]);
	K_ = K_*4^(m-2)
switch (mode)
case {'fw'}
	K_ref = [1 1 1; 1 -8 1; 1 1 1]/3;
case {'hw'}
	%K_ref = K;
	K_ref = [0,1,0; 1,-4,1; 0 1 0];
end
K_a = K_ - K_ref
K_a_ =K_a./2^(m-2)
sK_a_ = sum(K_a_)
K_./K_ref
pause

end

