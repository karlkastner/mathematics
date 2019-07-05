% 2016-03-07 16:15:26.221724969 +0100

%% skewness estimated by quantiles
function sk = qskewq(q,p)
	if (isvector(q))
		q = q';
	end
	sk = (q(1,:) - 2*q(2,:) + q(3,:))./(q(3,:)-q(1,:));
	% TODO, this constant is empirically determined
	% and slightly changes with the skewness
	%sk = 4.85*sk;
	%sk = 5.4763*sk - 5.8169*sk*sk;
% 	sk = 5.0733*sk - 19.1209*sk^3;
%	sk = 5.0695*sk - 18.9035*sk^3;
%
%	if (nargin()>1 && ~isempty(scaleflag) && scaleflag)
%		sk=4.9957*sk -15.1556*sk.^3;
%	end
end

