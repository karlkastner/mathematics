function plot_spectrum()

n=1000; E=[]; for idx=1:4; E = [E [1/idx^2 - 1./(idx+1:n).^2]]; end; E=1.09e7*E; idx=find(E>1e9/800); E=sort(E(idx)); lambda=1./(E(1:10))'*1e9;
clf
subplot(4,2,1)
for idx=1:length(lambda)
	l = lambda(idx);
	b = 380;
	r = 780;
	if (l >= 380 && l <= 780)
		c=4*[(l-b)/(r - b), 0 1-(l-b)/(r-b)].*(l - b)/(r - b);
%		line([l l],[0 1],'color', 1/64*c/max([1 1/64*c]),'Linewidth', 7 );
%		line([l l],[0 1],'color', 1/32*c/max([1 1/32*c]),'Linewidth', 6 );
%		line([l l],[0 1],'color', 0.6125*c/max([1 0.6125*c]),'Linewidth', 5 );
%		line([l l],[0 1],'color', 0.125*c/max([1 0.125*c]),'Linewidth', 4 );
%		line([l l],[0 1],'color', 0.25*c/max([1 0.25*c]),'Linewidth', 3 );
%		line([l l],[0 1],'color', 0.5*c/max([1 0.5*c]),'Linewidth', 2 );
%		line([l l],[0 1],'color', c/max([1 c]),'Linewidth', 1 );
		fadeline([l l],[0 1],c);
	end
end
set(gca,'FontName','serif')
xlabel('wavelength in nm', 'interpreter', 'latex');
%set( get( gca, 'XLabel' ), 'Interpreter', 'latex' ); 
%set(gca,'xtick',300:100:700)
%xtickabel(
set(gca,'Color',[0 0 0]);
set(gcf,'Color',[1 1 1]);
set(gca,'ytick',[]);
set(gcf, 'InvertHardCopy', 'off');
print -depsc hydrogen_spectrum_own.eps
%print -dpng hydrogen_spectrum_own.png

end

function fadeline(x,y,c)
	for idx=12:-1:0
		f = 2^-idx;
		p=1;
		c_ = f*c/max([1 f*c]);
		c_ = min(1,c_ + f*[1 1 1]*(1-norm(c_)/norm(c)));
		line(x,y,'color', c_,'Linewidth', idx+1 );
%		rectangle('Position', [x(1)-p*idx,y(1),2*p*(idx+1),y(2)-y(1)],'facecolor', c_, 'edgecolor',c_);
	end
end

