% Fri 24 Feb 13:33:38 CET 2023
function B = gbb_simulate(t,r,sigma,S0,varargin)
	S = gbm_simulate(t,r,sigma,S0,varargin{:});
	%B = S - (t-t(1))./(t(end)-t(1)).*(S(end,:)-S0);
	%B = exp(log(S) - (t-t(1))./(t(end)-t(1)).*log(S(end,:)./S0));
	%B = S.*exp(-(t-t(1))./(t(end)-t(1)).*log(S(end,:)./S0));
	t = cvec(t);
	B = S.*(S0./S(end,:)).^((t-t(1))./(t(end)-t(1)));
%	t  = t-t(1);
%	T  = t(end);
%	B_ = exp(log(S).*(1-t./T) - t./T.*log(S(end,:)./(S.*S0)));
%	s2 = exp(log(S).*(1-t./T) - t./T.*log(S(end,:)./(S.*S0)) - log(S0)).^2;
%	s2 = exp(log(S).*(1-t./T) - t./T.*log(S(end,:)./(S.*S0)) - log(S0)).^2;
%	rms(B-B_,'all')
end
