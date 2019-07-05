% 2013-08-16 14:14:10 UTC
% Karl Kästner, Berlin
%% convolution in 2d
function A_filt_= conv2_man(f1,f2,A,varagin)
	A_ = A;
	% 3d quick fix
	cl = class(A);
	for idx=1:size(A_,3)
	A = A_(:,:,idx);

	% padd with outermost cell
	A = [A(:,1)*ones(1,length(f2)-1,cl) A A(:,end)*ones(1,length(f2)-1,cl)];
	A = [ones(length(f1)-1,1,cl)*A(1,:); A; ones(length(f1)-1,1,cl)*A(end,:)];
%	A_filt = conv2(f1,f2,A,'valid');
	% filter
	A_filt = conv2(f1,f2,A,'same');
	% cut
	A_filt = A_filt(:,length(f2):end-length(f2)+1);
	A_filt = A_filt(length(f1):end-length(f1)+1,:);

	A_filt_(:,:,idx) = A_filt;
	end
end

