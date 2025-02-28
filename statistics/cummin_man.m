function [cm,mdx] = cummin_man(x)
	x = cvec(x);
	mdx = zeros(size(x));
	cm  = zeros(size(x));
	cm(1) = x(1);
	mdx(1) = 1;
	for idx=2:length(x)
		if (x(idx)<cm(idx-1))
			cm(idx) = x(idx);
			mdx(idx) = idx;
		else
			cm(idx)=cm(idx-1);
			mdx(idx)=mdx(idx-1);
		end
	end
end
