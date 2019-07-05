% 2016-05-03 20:00:35.525020683 +0200
%% transpose stacked 3x3 matrices
function At = transpose3(A)
	siz = size(A);
	if (length(siz)>2)
		siz3 = siz(3);
	else
		siz3 = 1;
	end
	At = zeros(siz(2),siz(1),siz3);
	for idx=1:siz(1)
	 for jdx=1:siz(2)
		At(jdx,idx,:) = A(idx,jdx,:);
	 end % for jdx
	end % for idx
end % transpose

