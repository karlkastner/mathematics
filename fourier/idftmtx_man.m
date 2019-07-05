% Thu  4 May 12:35:41 CEST 2017
% Karl Kastner, Berlin
%
%% inverse matrix for the discrete fourier transform in matlab style
%% with a limited number of columns, thus ignoring higher frequencies
%% keep 2nc+1 columns (mean and conj-complex pairs of nc frequencies)
%
function Fi = idftmtx_man(n,nc)
	nc = min(nc,floor(n/2));
	Fi = zeros(n,2*nc+1);
	Fi(1:n,1) = 1;
	id = (0:n-1)';
	I  = (1i).^(id*4/n);
	for idx=2:nc+1
		Fi(:,idx) = Fi(:,idx-1).*I;
		Fi(:,2*nc+3-idx) = conj(Fi(:,idx));
	end
	Fi=Fi/n;


end

