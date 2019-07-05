% Fri  5 Oct 15:38:10 CEST 2018
%
%% astar path finding alforithm
function [x,y,passed] = astar(val,x0,y0,p)
	% mouth to source: value must decrease
	x = [];
	y = [];
	passed = false(size(val));
	while (1)
		x=[x;x0]; 
		y=[y;y0];
		change = false; 
		x0_ = x0;
		y0_ = y0;
		passed(x0,y0) = true;
		v0 = -inf; % 0.5*val(x0,y0);
		%[x0,y0,v0]
		for dx=-1:1
		 for dy=-1:1
			% find next with maximum value
			if (~passed(x0_+dx,y0_+dy))
			if(val(x0_+dx,y0_+dy)>v0)
				% force descending
				if (val(x0_+dx,y0_+dy) <= val(x0_,y0_))
					x0 = x0_+dx;
					y0 = y0_+dy;
					v0 = val(x0,y0);
					change = true;
				end
			end % if
			end
		 end % dy
		end % dx
		if (~change)
			break;
		end
	end
end

