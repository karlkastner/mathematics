% Di 12. Jan 16:48:26 CET 2016
%% autocorrelation function
% c.f. genton 2000
function a = acf_genton(X,n)
	if (nargin() < 2)
		n = 20;
	end
	n = min(n,length(X)-1);
%	n = length(X)-1;
	a = zeros(n,1);
	for idx=1:n
	Qp2 = Q(X(1:end-idx+1) + X(idx:end))^2;
	Qm2 = Q(X(1:end-idx+1) - X(idx:end))^2;
	a(idx) =   (Qp2 - Qm2) ...
		./ (Qp2 + Qm2);
	end % for idx
end

function Q = Q(X)
	% TODO exact value
	k = 0.25;
	d = [];
	for idx=2:length(X)
		d = [d; abs(X(idx)-X(1:idx-1))];
	end
	Q = quantile(d,k);
end

