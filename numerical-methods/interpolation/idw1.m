% Do 14. Mai 10:04:38 CEST 2015
% Karl Kastner, Berlin
%
%% spatial average ny inverse distance weighting
function [val err] = idw2(X,Y,X0)
	% allocate memory
	val  = zeros(size(X0,1),1);
	err2 = zeros(size(X0,1),1);
	t    = tic();
    	for idx=1:length(X0)
%		if (toc(t) > 10)
%			disp(idx/length(X0));
%			t = tic();
%		end
    		d           = abs(X - X0(idx));
    		w           = 1./(d+1e-3);
    		w           = w/nansum(w.*isfinite(Y));
    		val(idx)    = nansum(w.*Y);
		dof         = sum(w).^2/sum(w.^2);
		err2(idx)   = nansum(w.*(Y - val(idx,1)).^2)/(dof-1);
    	end
	err = sqrt(err2);
end % idw2

