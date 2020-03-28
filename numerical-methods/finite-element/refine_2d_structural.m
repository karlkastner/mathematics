% Mon Jun 11 13:18:20 MSK 2012
% Karl Kästner, Berlin

function [X] = fem_2d_refine_structural(P,T,X,err,degen,order)
	X1 = X{1};
	X2 = X{2};
	n(1) = length(X1);
	n(2) = length(X2);
	dx = max([P(T(:,1),1) P(T(:,2),1) P(T(:,3),1)],[],2) - min([P(T(:,1),1) P(T(:,2),1) P(T(:,3),1)],[],2);
	dy = max([P(T(:,1),2) P(T(:,2),2) P(T(:,3),2)],[],2) - min([P(T(:,1),2) P(T(:,2),2) P(T(:,3),2)],[],2);

	% reshape the error
%	err1 = reshape(err.*(1),2*(n(2)-1),(n(1)-1));
%	err2 = reshape(err.*(1),2*(n(2)-1),(n(1)-1));
	err1 = reshape(err.*(dx>=dy),2*(n(2)-1),(n(1)-1));
	err2 = reshape(err.*(dx<=dy),2*(n(2)-1),(n(1)-1));
%	err1 = reshape(err.*(dx./dy).^2./degen,2*(n(2)-1),(n(1)-1));
%	err2 = reshape(err.*(dy./dx).^2./degen,2*(n(2)-1),(n(1)-1));
	% error inbetween four points is required, thus error of two triangles are added
	err1 = err1(1:n(2)-1,1:n(1)-1) + err1(n(2):2*(n(2)-1),1:n(1)-1);
	err2 = err2(1:n(2)-1,1:n(1)-1) + err2(n(2):2*(n(2)-1),1:n(1)-1);
	subplot(2,2,2)
	imagesc(err1)
	subplot(2,2,4)
	imagesc(err2)

	% max the error per row and column
	err2 = max(err2,[],2);
	err1 = max(err1,[],1)';
	%[XX YY] = meshgrid(X1(1:end-1),X2(1:end-2));
%	plot(X1(1:end-1), err1); hold on
%	plot(X2(1:end-1), err2,'r'); hold off
	drawnow
	pause(1)
	% find the maximum error
	err_max = max(max(err1),max(err2));

	% refine the mesh
	p = 0.5^order;

	% refine in x-direction
	X_ = zeros(2*n(1),1);
	n_ = 0;
	for idx=1:n(1)-1
		% take over old point
		n_ = n_+1;
		X_(n_,1) = X1(idx);
		if (err1(idx,1) > p*err_max)
			% add an additional centre point
			n_ = n_+1;
			X_(n_,1) = 0.5*(X1(idx,1) + X1(idx+1,1));
		end
	end % for idx
	% last point
	n_ = n_+1;
	X_(n_,1) = X1(n(1),1);
	% truncate the X_-vector to the correct length
	X{1} = X_(1:n_,1);

	% refine in y-direction
	X_ = zeros(2*n(2),1);
	n_ = 0;
	for idx=1:n(2)-1
		% take over the old point
		n_ = n_+1;
		X_(n_,1) = X2(idx);
		if (err2(idx,1) > p*err_max)
			% add an additional centre point
			n_ = n_+1;
			X_(n_,1) = 0.5*(X2(idx,1) + X2(idx+1,1));
		end % if
	end % for idx
	% last point
	n_ = n_+1;
	X_(n_,1) = X2(n(2),1);
	% truncate the X_-vector to the correct length
	X{2} = X_(1:n_,1);
end % fem_2d_refine_structural()

