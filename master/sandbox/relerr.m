function [nErr Err] = relerr(E)
	nErr = sqrt(sum( ((E-E(:,end)*ones(1,size(E,2)))./abs(E(:,end)*ones(1,size(E,2)))).^2 ));
	Err = sqrt(( ((E-E(:,end)*ones(1,size(E,2)))./abs(E(:,end)*ones(1,size(E,2)))).^2 ));
end

