% 2015-03-24 16:07:12.758711037 +0100
% Karl Kastner, Berlin
%
%% autoccorrelation of the columns of X 
function a = autocorr_man(x)
	if (isvector(x))
		x=cvec(x);
	end
	%for idx=1:size(X,2)
		x = bsxfun(@minus,x,mean(x));
		s20 = sqrt(mean(x(1:end-1,:).^2)*mean(x(2:end,:).^2));
		s21 = mean(x(1:end-1,:).*x(2:end,:));
		a = [ones(1,size(X,2)); (s21./s20).'];
	%end
end

