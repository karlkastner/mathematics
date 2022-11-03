% 2022-05-17 10:59:24.794758463 +0200
function fcut = fcut_gausswin(Lw)
	fcut = sqrt(-log(0.5)/4)./(pi*Lw);
end

