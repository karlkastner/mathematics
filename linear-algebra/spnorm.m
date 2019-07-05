% 2017-11-14 16:35:57.842755851 +0100
%% frobenius norm
function n = spnorm(A)
	n = sqrt(sum(A(:).^2));
end
