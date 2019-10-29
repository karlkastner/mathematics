% Oct 29  2010
% Karl Kästner, Berlin

% does not work due to matlab internal optimizaiton
% even with feature accel off
function s = sum_kahan(X)
%	s = zeros(size(X,1),1);
%	r = zeros(size(X,1),1);

%	for idx=1:size(X,2)
%		x = X(:,idx);
%		x = x + r;
%		s_new = s + x;
%		r = s_new - s - x
%		s = s_new;
%	end
	if (1 == isa(X,'single'))
		s = single(0.0);
		r = single(0.0);
	else
		s = 0.0;
		r = 0.0;
	end
	for idx=1:size(X,1)
		%x = X(:,idx) - r;
		x = X(idx);
		x = x - r;
		s_new = s + x;
		r = (s_new - s) - x;
		s = s_new;
	end	
end % function sum_kahan

