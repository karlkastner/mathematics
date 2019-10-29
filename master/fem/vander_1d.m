% Wed Jul 11 18:52:29 MSK 2012
% Karl Kästner, Berlin

function V = vander_1d(x,nv)
	V = zeros(size(x,1),nv);
	V(:,1) = 1;
	if (nv > 1)
		V(:,2) = x;
	end
	for idx=3:nv
		V(:,idx) = x(:,1).^(idx-1);
	end
end % vander_1d

