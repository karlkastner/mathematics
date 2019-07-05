% 2014-03-12 06:15:28.812847482 +0100
% Karl Kastner, Berlin
%
%% interpolate
function Vq = interp1(X,V,Xq,varargin)
	% set up the linear interpolation matrix
	jdx = 1;
	n  = length(X);
	nq = length(Xq);
%	buf = [[1:nq 1:nq]; ones(1,2*nq); NaN(1,2*nq)]';
	Vq = zeros(nq,size(V,2));
	% todo first nan values
	for idx=1:length(Xq)
		while(X(jdx+1) <= Xq(idx))
			jdx=jdx+1;
			if (jdx+1 > n)
				break;
			end
		end
			if (jdx+1 > n)
				% exact match
				if (Xq(idx) == X(jdx))
					buf(2*idx-1,:)     = [idx jdx 1];
				end
				break;
			end
		p = (X(jdx+1) - Xq(idx) )/(X(jdx+1)-X(jdx));
		buf(2*idx-1,:)     = [idx jdx   p];
		buf(2*idx,:)     = [idx jdx+1 1-p];
	end
	% sparse matrix multiplication
% sparse mm faster, but non sparse not
	%A = sparse(buf(:,1),buf(:,2),buf(:,3),nq,n);
	%Vq = A*V;
	for idx=1:size(buf,1)
		%Vq(idx,:) =   buf(2*idx-1,3)*V(buf(2*idx-1,2),:) ...
	        %            + buf(2*idx,  3)*V(buf(  2*idx,2),:);
		Vq(buf(idx,1),:) = Vq(buf(idx,1),:) +  buf(idx,3)*V(buf(idx,2),:);
	end
end
