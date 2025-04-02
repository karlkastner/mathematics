function s2 = powerpdf_var(a)
	s2 = (a-1)./((a-3).*(a-2).^2);
	s2(a<=3)=NaN;
end

