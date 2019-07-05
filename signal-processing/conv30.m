% 2014-07-14 10:25:29.657864377 +0200
% Karl Kastner, Berlin
%% convolve with rectangular window of lenght n
%% circular boundaries
% TODO rename this into filter
% TODO make n variable
% TODO 30 is even, better use odd number
function d_ = conv30(x)
	n = 15;
	if (1==ndims(x))
		d_ = [x(end-n+1:end); x(:); x(1:n)]; 
		d_ = conv(d_,ones(2*n,1)/(2*n),'same');
		d_ = d_(n+1:end-n);
	else 
		d_ = [x(end-n+1:end,:); x; x(1:n,:)]; 
		for idx=1:size(d_,2)
			d_(:,idx) = conv(d_(:,idx),ones(2*n,1)/(2*n),'same');
		end
		d_ = d_(n+1:end-n,:);
	end
end
