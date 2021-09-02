% Sun 11 Jul 14:02:01 CEST 2021
function R = imrotmat(n,a_rad)

todo corners cannot be saved, circular mask?

	x = (1:n(1))';
	x = repmat(x,1,n(2));
	x0 = flat(x);
	y = (1:n(2));
	y = repmat(y,n(1),1);
	y0 = flat(y);
%	x = x-mean(x);
%	y = y-mean(y);

	R = rot2(a_rad);

	xy = [x(:),y(:)];
	m  = mean(xy);
	xy = xy - m;
	% rotate
	xy = (xy*R');
	% shift back
	xy = xy+m;
	xy = round(xy);
	x = xy(:,1);
	y = xy(:,2);
	%x = x+n(1)/2;
	%y = y+n(2)/2;
%	x = round(x);
%	y = round(y);
	fdx=x<1;
	x(fdx) = x(fdx)+n(1);
	fdx = x>n(1);
	x(fdx) = x(fdx)-n(1);
	fdy=y<1;
	y(fdy) = y(fdy)+n(1);
	fdy = y>n(1);
	y(fdy) = y(fdy)-n(1);
	nn = prod(n);
	buf = [x0 + (y0-1)*n(1), x + (y-1)*n(1),ones(nn,1)];
	R = sparse(buf(:,1),buf(:,2),buf(:,3),nn,nn);
end
