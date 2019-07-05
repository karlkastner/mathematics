% 2016-03-15 14:20:20.258211991 +0100
function tV = interp_sn2(sS,sN,sdx,sV,tS,tN,tdx,S_range,L_max,order)
	tV = NaN(size(tS));
	sdx_ = unique(sdx);
	aspect = 1/3;
	for idx=1:length(sdx_)
		sdxi = sdx_(idx);
		fdx = sdx == sdxi;
		interp = TriScatteredInterp(aspect*sS(fdx),sN(fdx),sV(fdx),'linear');
		fdx = tdx == sdxi;
		tV(fdx,:) = interp(aspect*tS(fdx),tN(fdx));
	end
end


