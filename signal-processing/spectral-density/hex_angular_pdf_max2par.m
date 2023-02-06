% Sun 29 Jan 13:15:01 CET 2023
function k = hex_angular_pdf_max2par(Stc)
	lk = zeros(size(Stc));
	for idx=1:numel(Stc)
		lk(idx) = fzero(@(k_) hex_angular_pdf_max(exp(k_)) - Stc(idx),0);
	end
	k = exp(lk);
end

