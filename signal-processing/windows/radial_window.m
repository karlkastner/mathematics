% 2021-12-22 13:30:25.176493233 +0100
% Karl Kästner, Berlin
%
%% radial filter window in the 2d-frequency domain
%
function w = radial_window(r,rmax,method)
	if (nargin()<3)
		method ='gauss';
	end
	switch (method)
	case{0} % rect
		w = (r<=0.75*rmax);
	case {1} % linear 
		w    = min(1,max(2*(1-abs(r)./rmax),0));
	case {2} % cos
		%w = (1/2)*(1 - cos(2*pi*r./rmax));
		w = (1/2)*(1 - cos(2*pi*r./rmax)).*(abs(r)./rmax<1);
		w(abs(r)/rmax<0.5) = 1;
	case {'gauss'}
		s = rmax/sqrt(-log(0.25));
		w = normpdf(r,0,s)/normpdf(0,0,s);
	end % switch
end

