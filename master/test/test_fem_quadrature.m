% Fri Jul 13 14:18:50 MSK 2012
% Karl Kästner, Berlin

function test_fem_quadrature()

% 1D

C1_nc = {
'int_1d_nc_2',
'int_1d_nc_3',
'int_1d_nc_4',
'int_1d_nc_5',
'int_1d_nc_6'
}

C1_gauss = {
'int_1d_gauss_1',
'int_1d_gauss_2',
'int_1d_gauss_3',
'int_1d_gauss_4',
'int_1d_gauss_5',
'int_1d_gauss_6',
}

% 2D
C2_gauss = {
'int_2d_gauss_1',
'int_2d_gauss_3',
'int_2d_gauss_6',
'int_2d_gauss_7',
'int_2d_gauss_9',
'int_2d_gauss_12',
'int_2d_gauss_13',
'int_2d_gauss_16',
'int_2d_gauss_25',
}

C2_nc = {
'int_2d_nc_3',
'int_2d_nc_6'
'int_2d_nc_10',
'int_2d_nc_15',
'int_2d_nc_21',
};

% 3D

C3_nc = {
'int_3d_nc_4',
'int_3d_nc_6'
'int_3d_nc_11',
'int_3d_nc_20',
}
C3_gauss = {
'int_3d_gauss_1',
'int_3d_gauss_4',
'int_3d_gauss_5',
'int_3d_gauss_11',
'int_3d_gauss_14',
'int_3d_gauss_15',
'int_3d_gauss_24',
'int_3d_gauss_45',
};

%	test(C1)
%	test(C2)
%	test(C3)
	N = 2.^(1:6);

	subplot(2,3,1)
	tq(C1_nc,N,1)
	subplot(2,3,4)
	tq(C1_gauss,N,1);

	subplot(2,3,2)
	tq(C2_nc,N,2);
	subplot(2,3,5)
	tq(C2_gauss,N,2);

	N = 2.^(1:4);
	subplot(2,3,3)
	tq(C3_nc,N,3);
	subplot(2,3,6)
	tq(C3_gauss,N,3);

end

function test(C)
	for idx=1:length(C)
		[w b] = feval(C{idx});
		C{idx}
		e = ([double(sum(vpa(w))) min(sum(b,2)) max(sum(b,2))] -1)
		if (max(abs(e)) > 1e-12)
			pause
		end
	end
end

function [I M] = tq(C,N,d)
	for idx=1:length(C)
		C{idx}
		[w b] = feval(C{idx});
		for ndx=1:length(N)
			n=N(ndx);
			switch(d)
				case{1}
				f = 1;
				fun = @(Q) sin(Q(:,2));
				[P T Bc] = mesh_1d_uniform(n, pi);
				case{2}
				f = 4;
				fun = @(Q) sin(Q(:,2)).*sin(Q(:,3));
				[P T Bc] = mesh_2d_uniform([n n],[pi pi]);
				case {3}
				f = 24;
				fun = @(Q) sin(Q(:,2)).*sin(Q(:,3)).*sin(Q(:,4));
				[P T Bc] = mesh_3d_uniform(n,pi);
			end
			M(ndx,idx) = size(P,1);
			I(ndx,idx) = integrate(P,T,fun,w,b);
		end
	end
	err = I-f
	loglog(N,abs(err))
	legend(C)
		I
		M
end

function I = integrate(P,T,func,w,b)
	I = 0;
	d = size(P,2);
	for tdx=1:size(T,1)
		A = [ones(d+1,1) P(T(tdx,1:d+1),:)];
		area = 0.5*abs(det(A));
		Q = b*A;
		F = feval(func,Q);
		I = I + sum(area*w.*F);
	end
end

