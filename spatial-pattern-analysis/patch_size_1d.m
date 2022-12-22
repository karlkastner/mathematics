% 2022-07-02 17:33:56.983320797 +0200
function c = patch_size(b)
	s = size(b);
	c = zeros(s);
	for idx=1:s(2)
		% cont
		k = 1;	
		for jdx=2:s(1)-1
			if (~isnan(b(jdx,idx)))
			if (b(jdx,idx))
				k = k+1;
			else
				if (k>0)
					c(k,idx) = c(k,idx)+1;
				end
				k = 0;
			end % if b
			end % if isnan
		end % for jdx
	end % for idx
end % function

