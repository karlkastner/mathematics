% 2016-07-09 13:50:28.997251715 +0200 /home/pia/phd/src/lib/wrapphase.m
function a=wrapphase(a,o)
	if (nargin() < 2)
		o = 0;
	end
	o = o-pi;
	a = mod(a-o,2*pi)+o;
%	flag = a<-pi;
%	a(flag) = 2*pi+a(flag);
%	flag = a>pi;
%	a(flag) = 2*pi-a(flag);
end

