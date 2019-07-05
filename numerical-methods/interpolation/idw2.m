% Do 14. Mai 10:04:38 CEST 2015
% Karl Kastner, Berlin
%
%% spatial average by inverse distance weighting
function [val err] = idw2(X,Y,X0)
	% allocate memory
	val = zeros(size(X0,1),1);
	err2 = zeros(size(X0,1),1);
	t = tic();
    	for idx=1:length(X0)
		if (toc(t) > 10)
			disp(idx/length(X0));
			t = tic();
		end
    		d           = hypot((X(:,1) - X0(idx,1)),(X(:,2) - X0(idx,2)));
    		w           = 1./sqrt(d);
    		w           = w/nansum(w.*isfinite(Y));
    		val(idx,1)  = nansum(w.*Y);
		err2(idx,1) = nansum(w.*(Y - val(idx,1)).^2);	
    	end
	err = sqrt(err2);
end % idw2

