% So 3. Jan 11:46:50 CET 2016
%% bisection routine
% TODO golden ratio routine
% TODO introduce abstol / reltol
function [f x h] = line_search(func,g,hl,hr,fl,fr,it)
	if (false)
		f0 = fl;
		h = hr;
		while (true)
			x = h*g;
			f = func(-x);
			if (f < f0)
				break;
			end
			h = 0.5*h;
		end % while
	else

	hc = 0.5*(hl+hr);
	fc = func(-hc*g);
	if (it > 0)
		% bisect range around minimum in direction of the gradient
		if (fl < fr)
			[f x h] = line_search(func, g, hl, hc, fl, fc, it-1);
		else
			[f x h] = line_search(func, g, hc, hr, fc, fr, it-1);
		end
	else
		% end of iteration
		[f mdx] = min([fl fc fr]);
		switch (mdx)
		case {1}
			h = hl;
		case {2}
			h = hc;
		case {3}
			h = hr;
		end
		x = h*g;
	end

	end
end % line_search

