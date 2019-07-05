% Sat  6 May 10:44:34 CEST 2017
%% inverse of the short time fourier transform
%
% TODO at the moment prediction is only possible at sample points
function [val, ct, vf, obj] = itransform(obj,time)

	% coefficients
	c = flat(obj.coeff);

	% set up regression matrix
	if (nargin()>1 && ~isempty(time))
		[A1 A2] = obj.stftmat(time);
			%obj.t0,obj.dt,obj.nt,obj.T,obj.Ti,obj.order);
	else
		% reexpand at sampled points
		A1 = obj.A1;
		A2 = obj.A2;
	end
	
	% predict
	ct  = A2*c;
	val = A1*ct;
	if (nargout()>1)
		ct  = reshape(ct,obj.ni,[]);
	end
	% for individual frequency components
	if (nargout()>2)
		% I  = speye(obj.ni);
		c_ = zeros(size(A1,2),length(obj.T)+1);
		c_(1:obj.ni:end,1) = ct(1,:).';
		for idx=2:length(obj.T)+1
			c_(2*idx-2:obj.ni:end,idx) = ct(2*idx-2,:).';
			c_(2*idx-1:obj.ni:end,idx) = ct(2*idx-1,:).';
		end
		%vf = A1*(kron(ct,I))';
		vf = A1*c_;
	end

	% time
	%time  = obj.t0 + obj.dt*(0:obj.nt-1)';

	if (0)

	Ti = obj.Ti;
	T  = obj.T;
	t0 = obj.t0;

	% samples per interval
	n = round(Ti/dt);

	% number of intervals
	m = size(obj.coeff,2);


	A = fouriermtx(timei,T);

	c = obj.coeff;
	val = zeros(n*m,size(c,3));
	val_ = [];
	for cdx=1:size(c,3)
		vali = A*c(:,:,cdx);
		val_(:,:,cdx) = vali;
		val(:,cdx) = vali(:);
	end

	end
end % itransform

