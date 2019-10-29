function errmat = errmat(E, e)
	errmat = zeros(size(E));
	for idx=1:size(E,2)
		for jdx=1:length(e)
			errmat(jdx,idx) = min(abs( (E(:,idx) - e(jdx))/e(jdx) ) );
		end
	end
end

