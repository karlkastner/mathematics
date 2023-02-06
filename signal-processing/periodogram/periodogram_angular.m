% 2021-06-21 22:52:10.758702394 +0200
function [St,theta,A] = periodogram_angular2(Shat,L)
	n     = size(Shat);
	if (abs(L(1)-L(2))>sqrt(L(1)*L(2)*eps) || n(2) ~= n(1))
		error('matrix must be square')
	end

	% discretizing the maximum circumference ensures that sectors
	% are only 1 pixel wide
	m     = round(pi*n(1));
	
	theta = 2*pi*(0:m-1)'/m-pi;
	nr    = floor(n(1)/2);
	r     = (1:nr)';
	s     = sin(theta);
	c     = cos(theta);
	ind_Sr = [];
	buf = [];
	% the low frequency bins are part of several sectors and thus are weighted down in proportion
	w   = r./(pi*sum(r));
	x = cvec(r)*rvec(c);
	y = cvec(r)*rvec(s);
	% bilinear interpolation
	% note, this really has to be floor not fix
	i = floor(x);
	p = 1-(x-i);
	j = floor(y);
	q = 1-(y-j);

	buf = [  i(:),j(:)  ,flat(p.*q.*w),
                i(:)+1,j(:)  ,flat((1-p).*q.*w),
		i(:)  ,j(:)+1,flat(p.*(1-q).*w),
                i(:)+1,j(:)+1,flat((1-p).*(1-q).*w)];
	ind_Sr = repmat(flat(repmat((1:m),nr,1)),4,1);
	% this translates indices in the same mannler like ifftshift
	buf(:,1) = fourier_freq2ind(buf(:,1),n(1));
	buf(:,2) = fourier_freq2ind(buf(:,2),n(2));
	% indices in flattened matrix
	ind_S2d = sub2ind(n,buf(:,1),buf(:,2));	
	% averaging matrix
	A = sparse(ind_Sr,ind_S2d,buf(:,3),m,n(1)*n(2));
	% average spectrum
	St = A*flat(Shat);
	dt = 2*pi/m;
	St = St./(sum(St)*dt);
end

