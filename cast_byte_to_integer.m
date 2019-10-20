% Sun Jul 13 12:16:49 WIB 2014
%% cast byte to integer
function int = cast_byte_to_integer(b)
		if (isvector(b))
			b = rvec(b);
		end
	         int =        1*uint32(b(:,1)) ...
                     +      256*uint32(b(:,2)) ...
                     +    65536*uint32(b(:,3)) ...
                     + 16777216*uint32(b(:,4));
end

