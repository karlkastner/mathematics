if (0)
clf
subplot(2,2,1)
imagesc((abs(fft2(mask)).^2))
%mean(Shat,3))
colorbar
subplot(2,2,2)
imagesc(w_);
%imagesc(mean(Sbar,3))
colorbar
subplot(2,2,3)
imagesc(mean(Shat./Sbar,3))
colorbar
subplot(2,2,4)
imagesc(mean(Shat,3)./mean(Sbar,3))
colorbar
end
if (0)
	q = sort(ratio(:));
	n = length(q);
	p = (1:n)'/(n+1);

	mu = mean(q);
	sd = std(q);
	[d1,d2] = fisher_moment2param(mu,sd);
	
	Shat_  = zeros(size(mask));
	Shat__  = abs(fft2(mask)).^2;
	id = 1;
	Shat_(1,1) = Shat__(1,1);
	Shat_bar = circfilt2(Shat_,nf);
	ratio_max = Shat_(1,1)./Shat_bar(1,1);
clf
imagesc(Shat_)
Shat_(1,1)
Shat_bar(1,1)
max(ratio(:))./ratio_max
pause

% nb: it is not quite clear what the generalization should look like, as ratio has
%     to be devided by nf^2 to match a beta distribution
%	[a,b] = beta_moment2param(mu,sd)

	stat.f.d1 = d1;
	stat.f.d2 = d2;
	stat.mu = mu;
	stat.sd = sd;
end

	% determine the maximum value the distribution can take
	


if (0)
mu
%beta_mean(a,b)
sd
%beta_std(a,b)
%q_  = betainv(p,a,b)
q__ = finv(p,d1,d2)
figure(1)
clf
%plot(p,[q,q__,q_])
length(p)
d_ = round(length(p)/100);
plot(p(1:d_:end),[q__(1:d_:end)])
drawnow
'honk'
pause
end


