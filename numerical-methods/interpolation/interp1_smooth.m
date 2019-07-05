% 2014-08-24 12:32:45 +0200
% Karl Kastner, Berlin
%
function [Vi dVi ddVi] = interp1_smooth(X,V,Xi,r0,order)
	% polynomial derivative scale factors
%	one = ones(1,size(V,2));
%	d = ((0:order)'*one);
%	dd = d.*((-1:order-1)'*one);
	dc  = zeros(order+1,size(V,2));
	ddc = zeros(order+1,size(V,2));
	for idx=1:length(Xi)
		% TODO, speed up by qtree
		D2 = (X - Xi(idx)).^2;
		%[mv mdx] = min(D2);
		fdx = find(D2 < r0*r0);
		Xs = X(fdx);
		mu = mean(Xs);
		Xs = Xs-mu;
		A = [];
		Ai = [];
		%W = (sqrt(D2(fdx)));
		W = cos(0.5*pi*sqrt(D2(fdx))/r0);
		W = W/sum(W);
		W = diag(W);
		% TODO preallocate A and Ai
		for jdx=1:order+1
			A  =  [A Xs.^(jdx-1)];
			Ai  = [Ai (Xi(idx)-mu).^(jdx-1)];
		end
%		c = A \ V(fdx,:);
		% TODO, use better chol here or indicate at least symmetry to the solver
		c = (A'*W*A) \ (A'*W*V(fdx,:));
		Vi(idx,:) = Ai*c;
		if (nargout() > 1)
			%for jdx=1:order
			%	dAi(:,jdx) = jdx*Ai(:,jdx);
			%end
			%dc = c.*d;
			dc(1,:) =   c(2,:);
			dc(2,:) = 2*c(3,:);
			dc(3,:) = 3*c(4,:);
			dVi(idx,:) = Ai*dc; %(2:end,:);
		end
		if (nargout() > 2)
			%for jdx=1:order-1
			%	ddAi(:,jdx) = jdx*(jdx+1)*Ai(:,jdx);
			%end
			%ddVi(idx,:) = ddAi*c(3:end,:);
			%ddc = c.*dd;
			ddc(1,:) = 2*c(3,:);
			ddc(2,:) = 6*c(4,:);
			ddVi(idx,:) = Ai*ddc;
		end
	end % for idx
end % interp1_smooth()

