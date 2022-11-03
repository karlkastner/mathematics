% Mon 13 Jun 13:50:39 CEST 2022
% gelbaum 2014
% c.f wang 2008, p 48 -> the definition of wang is different (wrong?)
%
%
% dietrich, newsham: only for symmetric covariances, which is not the case
function y = brownian_motion_2d_fourier(n,kmax)
	if (length(n)<3)
		n(3)=1;
	end
	if (nargin()<2)
		kmax=inf;
	end
	e = randn(n(1)*n(2),n(3));
	y = zeros(n(1)*n(2),n(3));
	x1 = innerspace(0,1,n(1))';
	x2 = innerspace(0,1,n(2))';
	k = 1;
	s = 1.5;
	for idx=1:min(n(1),kmax)
			v1 = sqrt(2)*sin((idx-1/2)*pi*x1);
			l1 = (idx-1/2)^2*pi^2;
		for jdx=1:min(n(2),kmax)
			v2 = sqrt(2)*sin((jdx-1/2)*pi*x2);
			l2 = (jdx-1/2)^2*pi^2;
			sl = (l1 + l2).^(s/2);
			%l = l1*l2;
			v = kron(v1,v2); %/sqrt(n(1)*n(2));
			y = y + 1/sl*v*e(k,:);
			k = k+1;
		end
	end
	y = y; %*sqrt(n(1)*n(2));
	y = reshape(y,n);
end
