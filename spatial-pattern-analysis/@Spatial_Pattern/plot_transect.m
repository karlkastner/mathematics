% Tue 20 Jul 23:47:34 CEST 2021
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% plot 1D pattern
%
function transect_plot(obj,meta)
	% plots
	fflag = meta.pflag;

	b    = obj.b;
	R    = obj.R;
	S    = obj.S;
	stat = obj.stat;

	printf('Periodogram Statistics\n');
	printf('max Sb %f\n',stat.Sc.bartlett); %*fx(mdx));
	printf('f_mean %f\n',stat.f.mean);
	printf('f_std %f\n',stat.f.std);
	printf('sd_k %f\n',stat.f.std);

	% workaround for matlab bug mistaking pi as variable name
	pi = 3.146;
	printf('k_c BM: %g\n',2*pi*stat.fc.brownian);
	printf('lambda_c : %g\n',stat.wavelength.brownian)
	printf('Sc : %g\n',stat.Sc.brownian);
	printf('periodicity test p-value%f\n',stat.p_periodic);
%	printf('p inclusive %f\n',pi);

	printf('stationarity test p : %f\n', stat.p_stationary);
	printf('stationarity test D : %f\n', stat.D_stationary);

	xscale = stat.fc.ref;
	if (isfield(meta,'scale'))
	switch (meta.scale)
	case {'L'}
		yscale = 1/obj.L;
	case {'lambda'}
		yscale = xscale;
	otherwise
		error('here');
	end
	else
		yscale = xscale;
	end
	% pattern
	x   = obj.x;
	[fx,fdx]  = obj.fx; 
	%fdx = fx>=0;
	splitfigure([2,2],[1,1],fflag);
	cla();
	drawnow();
	fdx_ = xscale*(x-x(1))<meta.pattern.xlim(2);
	b_ = b(fdx_);
	b_ = b_-min(b_);
	b_ = b_/rms(b_);
	plot(xscale*x(fdx_),b_,'linewidth',1.5)
	xlim(meta.pattern.xlim);
	ylim([0,3.5]);
	set(gca,'ytick',meta.pattern.ytick);
	xlabel(meta.pattern.xlabel,'interpreter','latex');
	ylabel(meta.pattern.ylabel,'rot',90,'interpreter','latex');
	drawnow();

	% periodogram and density
	splitfigure([2,2],[1,2],fflag);
	cla();
	drawnow();
	if (stat.p_periodic < obj.opt.confidence_level)
	stem(fx(fdx)/xscale,yscale*S.hat(fdx),'-','linewidth',4,'marker','none');
	ylim([0, 1.05*max(yscale*S.hat(fdx))])
	else
	errorarea2(fx(fdx)/xscale,yscale*[S.c(fdx,1),zeros(size(S.c(fdx),1),1),S.c(fdx,2)],meta.areacol,'FaceAlpha',1)
	hold on;
	plot(fx(fdx)/xscale,yscale*S.hat(fdx),'k-','linewidth',1);
	plot(fx(fdx)/xscale,yscale*S.ref(fdx),'r','linewidth',1.5);
	ylim(meta.periodogram.ylim*yscale/xscale);
	end
	%set(gca,'ytick',meta.periodogram.ytick);
	xlabel(meta.periodogram.xlabel,'interpreter','latex')
	ylabel(meta.periodogram.ylabel,'rot',0,'interpreter','latex');
	xlim(meta.periodogram.xlim)
	vline((1:10),'linestyle','--','color','k')
	drawnow();

	% autocorrelation
	splitfigure([2,2],[1,3],fflag);
	cla();
	drawnow();
	if (stat.p_periodic < obj.opt.confidence_level)
	R.rate = acf_decay_rate(x,R.hat);
	plot(x*xscale, R.hat, 'linewidth',1.5);
	hold on;
	h = plot(x*xscale,exp(-R.rate*x),'k--','linewidth',1.5);
else
	plot(x*xscale,R.ref,'linewidth',1.5);
	hold on
	if (~isempty(stat.acf_decay_rate))
	h   = plot(x*xscale,exp(-stat.acf_decay_rate*x),'k--','linewidth',1.5);
	end
end
	h.HandleVisibility = 'off';
	xlim(meta.acf.xlim);
	ylim(meta.acf.ylim);
	ytick(meta.acf.ytick);
	xlabel(meta.acf.xlabel,'interpreter','latex');
	ylabel(meta.acf.ylabel,'rot',0,'interpreter','latex');
	hline(0);
	drawnow()
	

	% QQ-plot
	fdx = (fx>=0);
	m = sum(fdx);
	pq = (1:m)'/(m+1);
	a = 1;                                                          
	b_ = (stat.m-1);                                                     
	qbeta = stat.m*betainv(pq,a,b_);
	d1 = 2;
	d2 = 2*(stat.m-1);
	qfisher = finv(pq,d1,d2);
	qchi2 = 1/2*chi2inv(pq,2);

	% note : for bartlett, the distribution is not beta,
	%        only close to beta for a small number of splits
	% the distribution is exactly beta for smoothing
	if (true) %nargin()<6)
		% empirical
		qref = qbeta;
		S.ref = S.bartlett;
		S.ref = S.filt;
		%S.ref = qfisher;
	else
		% exact
		qref = qchi2;
		S.ref = S.ref;
	end

	%figure(200)
	%clf();
	splitfigure([2,2],[2,1],fflag);
	cla();
	plot(pq,1/2*chi2inv(pq,2));
	hold on;
	plot(pq,sort(S.hat(fdx)./S.flat(fdx)));
	plot(pq,sort(S.hat(fdx)./S.bartlett(fdx)));
if (0)
	plot(pq,sort(S.hat(fdx)./S.mui(fdx)));
end
	
	splitfigure([2,2],[2,2],fflag);
	cla();
	h = plot([0,20],[0,20],'linewidth',0.5);
	h.HandleVisibility='off';
	hold on;
	plot(qchi2,sort(S.hat(fdx)./S.flat(fdx)),'.k','linewidth',1.5);
	leg_C = {'flat'};
	%plot(qbeta,sort(S.hat(fdx)./S.mui(fdx)),'.b','linewidth',1.5);
	plot(qref,sort(S.hat(fdx)./S.ref(fdx)),'.r','linewidth',1.5);
	leg_C{end+1} = 'empirical';
	if (0)
	plot(1/2*chi2inv(p,2),sort(S.hat(fdx)./S.f(fdx,2)),'.b','linewidth',1.5);
	leg_C{end+1} = {'lorentzian'};
	end
	xlabel('expected quantile');
	ylabel('sample quantile');
	axis(8*[0,1,0,1]);
	axis equal;
	axis square;
	legend(leg_C{:},'location','southeast');
	set(gca,'colororder',meta.colororder);
	
	% stationarity
	%figure(300);
	%clf
	nb=5;
	splitfigure([2,2],[2,3],fflag);
	cla
	try
		q = periodogram_qq(b,nb);
		plot(q(:,1),q(:,2),'.');
		hold on;
		mq=1.05*(max(q(:)));
		plot([0,mq],[0,mq],'k');
		axis square;
		axis([0,mq,0,mq])
		xlabel('expexcted')
		ylabel('sampled');
	catch e
		e
	end
	
	% average band shape
	splitfigure([2,2],[2,4],fflag);
	cla
	nf = 100;
	[y_mu,y_sd,yy] = average_wave_shape(b,nf);
	y_mu = circshift(y_mu,+15);
	y_mu = trifilt1([y_mu;y_mu;y_mu],21);
	y_mu = y_mu(nf+1:2*nf);
	y_mu = y_mu/rms(y_mu);
	n = length(y_mu);
	x = linspace(0,1,nf)';
	plot(x,y_mu,'linewidth',1.5);
	ylim([0,max(y_mu)*1.02])
	%xlabel(meta.pattern.xlabel);
	%ylabel(meta.pattern.ylabel);
	%plot(([circshift(y_sd./sqrt(size(yy,2)),+11),circshift(y_mu,+11)]./mean(y_mu)))
	
	% QQ
	if (meta.pflag)
		ps = meta.plotscale;
		aspect = meta.aspect;
		basename = obj.opt.basename;
		pdfprint(11,[basename,'-pattern.pdf'],ps,aspect);
		pdfprint(12,[basename,'-psectral-density.pdf'],ps,aspect);
		pdfprint(13,[basename,'-autocorrelation.pdf'],ps,aspect);
		pdfprint(22,[basename,'-qq.pdf'],ps,aspect);
%	pdfprint(100,'img/observed-pattern-stationarity-qq.pdf',sp);
	%	pdfprint(24,[basename,'-shape.pdf'],sp);
	end

