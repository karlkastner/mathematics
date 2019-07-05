% Wed Oct 29 11:50:46 CET 2014
% Karl Kastner, Berlin
%
%% jacknife estimate of mean, bias and standard error
%% theta0 : estimate from all samples
%% thetad : set of estimates obtained by leaving out one data point each
%%          last dimension of theta is assumed to be the jackknife dimension
function [mu, bias, serr2, C, N] = estimated_STATIC(theta0,thetad,d,hat)
	s    = size(thetad);
	% fix for 1D
	if ( 2 == length(s) && 1 == s(2) )
		jdim = 1;
	else
		jdim = length(s);
	end

	% n the number of jacknife pseudo samples (only true for d = 1)
	n    = size(thetad,jdim);
	N    = sum(isfinite(thetad),jdim);
	% according to Simonoff 1986, the JN has to be reweighed in case
	% of onlinear functions, as otherwise the bias get's worse
	mt   = nanmean(thetad,jdim);
	if (nargin() < 4 || isempty(hat))
		mu   = (N/d).*theta0 - ((N-d)/d).*mt;
	else
		% weigh down of extreme influence points (RLQ Simonoff 1986)
		% this is only for the bias, not for the variance
		w  = rvec((1-diag(hat)));
                %w = w/sum(w);
		%mt = nansum(repmat(w,[s(1) 1]).*thetad,jdim);
		%mt = nansum(repmat(w,[1 s(2)]).*thetad,jdim);
		%mt = nansum(repmat(w,[1 s(2) s(3)]).*thetad,jdim);
		n_ = length(w);
		for idx=1:length(w)
			P(:,idx) = theta0 - n_*w(idx).^2*(theta0 - thetad(:,idx));
		end
		mu = squeeze(mean(P,2)); %jdim-1));
	end
	bias = theta0 - mu;

	% fix for singleton dimension
	if (1 == jdim)
		s = [n, 1];
	else
		s = [ones(1,jdim-1) n];
	end
	% error variance - diagonal of the covariance matrix
	if (isempty(hat))
		serr2 = (N-d).*nanmean((thetad - repmat(mt, s)).^2,jdim);
	else
		n = size(thetad,jdim);
		% construct pseudo variables
		for idx=1:length(w)
			P(:,idx) = theta0 + n*w(idx)*(theta0 - thetad(:,idx));	
		end
		% mean of pseudo variables
		thetaj = mean(P,2);
		% deviation of pseudo variables from their mean
		D = squeeze(P - repmat(thetaj,1,n));
		% variance
		serr2 = 1/(n*(n-1))*sum(D.^2,2);
	end
	% covariance matrix
	if (nargout() > 3)
		% pseudo variables
		n = size(thetad,jdim);
		% for the variance, the samples are weighed with w, not with w^2
		% according to simonoff
		if (nargin() < 4 || isempty(hat))
			%P = squeeze(n*repmat(theta0,s) - (n-1)*thetad);
			P = (n*repmat(theta0,s) - (n-1)*thetad);
			if (size(P,2) == 1)
				P = reshape(P,size(P,1),size(P,3));
			end
			%D = P - repmat(mean(P,2),1,n);
		else
			for idx=1:length(w)
				P(:,idx) = theta0 + n*w(idx)*(theta0 - thetad(:,idx));
			end
		end
		% mean of pseudo variables
		thetaj = mean(P,2); %jdim);
		% deviation of pseudo variables
		D = squeeze(P - repmat(thetaj,1,n)); % s
		% covariance matrix
		C = zeros(length(theta0));
		for idx=1:n
			C = C + D(:,idx)*D(:,idx)';
			%C_ = C_ + D_(:,idx)*D_(:,idx)';
		end
		%C_ = 1/n*C_;
		% C_ seems factor 3 to small
		C = 1/(n*(n-1))*C;
	end
	
%	n    = size(theta,jdim);
%	mt   = mean(theta,jdim);
%	mu   = n/d*theta0 - (n-d)/d*mt;
%	verr = (n-d)*mean((theta - repmat(mt,[ones(1,jdim-1) n])).^2,jdim);
%	serr = sqrt(mean(res2));
end % jacknife_from_data

