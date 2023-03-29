% Fri  4 Mar 14:19:02 CET 2022
function c = cwt_man(y,L0,dx,w)
	for idx=1:length(L0)
		xf  = linspace(-0.5*w*L0(idx),0.5*w*L0(idx),round(w*L0(idx)/dx))';
		win = triwin(xf);
		%win = rectwin(xf);
		%win = hanwin(xf);
		%win = ones(length(xf));
		win = win/sum(win);
		fc  = 2*win.*cos(2*pi*xf/L0(idx));
		fs  = 2*win.*sin(2*pi*xf/L0(idx));
		f   = fc+1i*fs;
		c(:,idx) = conv(y,f,'same');
	end
end

