% Mi 14. Okt 15:48:04 CEST 2015
% Karl Kastner, Berlin
%
%% non-central higher order moments of the multivariate normal distribution
%%
%% c.f. Moments and cumulants of the multivariate real and complex Gaussian distributions
%% 
%% note : there seem to be some typos in the original paper, 
%%	for x^4 cii^2, the square seems to be missing
%% mu : nx1 mean vector
%% C  : nxn covariance matrix
%% k  : nx1 powers of variables in moments
function M = comoment(mu,C,k)
	% expand k
	K = [];
	for idx=1:length(k)
	for jdx=1:k(idx)
		K = [K;idx];
	end
	end
	n = length(K);
	M = 0;
	% choose first (fixed to 1)
	% choose second
	for idx=2:n
		% choose third
		for jdx=2:n
			if (jdx~=idx)
			% choose fourth
			for kdx=2:n
				if (kdx~=idx && kdx ~= jdx)
%				[1 idx jdx kdx]
				M = M + C(K(1),K(idx))*C(K(jdx),K(kdx));
			end
			end % kdx
			end
		end % jdx
	end % idx
%	M_ = choose(K,[],length(K));
%	M = 0;
%	for idx=1:size(M_,1)
%		m = 1;
%		for jdx=1:2:size(M_,2)
%			% TODO, for 6 this become triple, for two this is just 1 factor
%			m = m*C(K(M_(idx,jdx)),K(M_(M_(idx,jdx+1))));
%		end
%		M = M + m;
%	end
end % function

function M = choose(K,banned,level)
	if (level > 0)
		M = [];
		for idx=1:length(K)
			if (~ismember(idx,banned))
				M = [M; choose(K,[banned idx],level-1)];
			end
		end
	else
		M = banned;
	end
end

