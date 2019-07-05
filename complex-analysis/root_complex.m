% 2017-12-04 11:54:45.139507162 +0100
%% root of a complex number
function [a, b] = root_complex(re,im)
	% this fails if a==0,which is the case if im == 0 and re < 0
	a=sqrt(1/2*(re + sqrt(im.^2+re.^2)));
	b=im./(2*a);
	if (issym(re) || issym(im))
		fdx = isAlways(im == 0) & isAlways(re < 0);
	else
		fdx = (im == 0 & re < 0);
	end
	b(fdx) = sqrt(-re(fdx));
end

