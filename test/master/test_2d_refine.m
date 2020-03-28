function test_2d_refine(disp)

path(path,'fem')

clf
[P T B] = mesh_2d_uniform(3,1);
[T n] = restore_cw(P, T, B);
check_area_2d(P,T,B)
G = zeros(size(T,1),1);

%T = [1 2 3];
%B= [1 2 1; 2 3 2; 3 1 3];
[Nm ] = fem_2d_element_boundary(P, T, B);

%test with negative/mixed coords!
k=10;
for idx=1:k
	l=size(T,1);
	'size'
	[size(P,1) size(T,1) size(B,1)]
	ex=0;
	% check the area
	'check area'
	[c l_ d] = check_area_2d(P,T,B);
	if (abs(c-1) > 1e-12 || abs(l_-4) > 1e-12 || d > 0)
		'area failed'
		ex = 1;
	end
	'restore cw'
	[T_ n_] = restore_cw(P, T, B);
	if (n_ > 0)
		'cw failed'
		ex=2;
	end

	'neighbourhood test'
	try 
	        [Nm_ ] = fem_2d_element_boundary(P, T, B);
		'pass'
		[idx_ jdx] = find(Nm-Nm_ ~= 0);
		N = [Nm Nm_ Nm-Nm_];
		if (norm(Nm - Nm_) > 0)
			'neighboorhood failed'
			ex=3;
			[unique(idx_) N(unique(idx_),:)]
		end
	catch err
		'neighboorhood exception'
		ex = 4;
		N = [];
		idx_ = [];
	end
if (nargin > 0 && 0 ~= disp)
	figure(idx+1); clf;
end
	if (ex > 0)
	%	[(1:length(T))' T]
		subplot(1,1,1); clf
		display_2d(P,T,B,3); axis equal; axis tight
		hold on
		plot( mean(P(T(idx_,:),1),2),mean(P(T(idx_,:),2),2) , 'r.');
		%diff(sort(10*P(:,1) + P(:,2)))'
		drawnow
		error('test','here')
	end


%	display_2d(P,T,B,3); axis equal; axis tight
		
	l=size(T,1);
	M=unique(randi(l,max(1,round(l/16)),1));
	[P T B Nm G] = fem_2d_refine(P,T,B,M,Nm,G);
%	[P T B Nm] = flip(1:size(T,1),P,T,B,Nm);

	if (nargin > 0 && 0 ~= disp)
		subplot(2,2,4);
		display_2d(P,T,B,7); axis equal; axis tight
		title('after flip')
	end

%	T
%	display_2d(P,T,B);
%	[c l d] = check_area_2d(P,T,B);
%	if (abs(c-4) > 1e-12 || abs(l-8) > 1e-12 || d > 0)
%		'error'
%		cla
%		c
%		l
%		d
%		axis([-0.1 2.1 -0.1 2.1])
%		display_2d(P,T,B,(1:l)'/l,2);
%	else
%		'ok'
%	end
%	subplot(2,k,idx);
%display_2d(P,T,B) %,0.5*ones(l,1));
%	axis([-0.1 2.1 -0.1 2.1])

end

end




