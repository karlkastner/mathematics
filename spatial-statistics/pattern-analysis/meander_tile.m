function [z1] = meander_tile(n1)

	% prototype tile
	z1 = zeros(n1,n1);
	dir = [ 0,1
                1,0
                0,-1
               -1,0];
                
	z1(1,1:n1-1) = 1;
	id = [1,n1-1];
	ddx = 1;
	n1_ = n1(1)-2;
	while (n1_>0)
		ddx = mod(ddx,4)+1;
		id_(1) = id(1)+n1_*dir(ddx,1);
		id_(2) = id(2)+n1_*dir(ddx,2);
		if (dir(ddx,1)~=0)
			z1(id(1):dir(ddx,1):id_(1),id(2)) = 1;
		else
			z1(id(1),id(2):dir(ddx,2):id_(2)) = 1;
		end
		id = id_;
		n1_ = n1_ - 1;
	end
end

