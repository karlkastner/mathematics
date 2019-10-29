% Fri Sep 23 13:31:35 MSD 2011

function C = mmmul_accurately(A,B)
	% test sizes
	if (size(A,2) ~= size(B,1))
		'error, matrix sizes do not match'
	end
	% allocate memory
	C = zeros(size(A,1), size(B,2));
	% multiply
	for idx=1:size(B,2)
		for jdx=1:size(A,1)
			%C(jdx,idx) = A(jdx,:)*B(:,idx);
			C(jdx,idx) = pairwise_sum(A(jdx,:).*B(:,idx)');
		end
	end
end % mmul_accurately


%factorial
% n=factorial(6)+1; tic; [V E] = eig(full(hydrogen_boxed(n,1))); toc, np=1; for idx=1:np; subplot(sqrt(np),sqrt(np),idx); [s sdx] = sort(diag(E)); plot(V(:,sdx(end+1-idx))); ylim([-0.15 0.15]); xlim([1 n]); end


function test_sum(k)
	for idx=1:k
		n=2^idx;
		r=rand(n,1).*rand(n,1);
		rs = single(r);
		r_ = sum(r)
		%[n sum(rs)-r_ l_sum(rs)-r_ f_sum(rs)-r_ v_sum(rs)-r_]
		[n sum(rs)-r_ l_sum(rs)-r_ v_sum(rs)-r_]
	end
end



function x = f_sum(x)
	l=length(x);
	s = 0;
	for idx=1:l
			s = s + x(idx);
	end
	x = s;
	%x = ones(1,l)*x;
end

function x = v_sum(x)
	l=length(x);
%	s = 0;
%	for idx=1:l
%			s = s + x(idx);
%	end
%	x = s;
	x = ones(1,l)*x;
end
