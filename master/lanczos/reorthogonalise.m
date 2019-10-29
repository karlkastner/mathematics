function Q = reorthogonalise(Q)
	for idx=2:size(Q,2)
		Q(:,idx)=mgs(Q,[],[],Q(:,idx),idx-1);
	end
end

