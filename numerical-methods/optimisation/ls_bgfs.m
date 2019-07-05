% Tue 16 Aug 15:05:53 CEST 2016
% Karl Kastner, Berlin
%% least squares by the bgfs method
function lsbgfs()
	B = eye(p);
	while (1)
		[f g] = afun();
		s = -(B\g);
		% step
		% TODO, line search
		x = x + s;

		case {'bgfs'}
			s = -(B\g);
			% step
			% TODO, line search
			x = x + s;
	
			dx = x - xold;
			dg = r - rold;
			% BGFS update
			B = bgfs(B,dx,dg);
		end


		dx = x - xold;
		dg = g - gold;
		% BGFS update
		B = bgfs(B,dx,dg);

	end
end

function B = bgfs(B,ds,dg)
	Bdx = B*dx;
	B = B + (1/(dx'*dg)*dg)*dg' - (1/(dx'*Bdx))*(Bdx*Bdx');
end

