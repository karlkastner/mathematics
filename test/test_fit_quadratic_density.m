% Thu  6 Oct 14:38:42 CEST 2022

f  = linspace(0,1e2,1e2);
fc = 5;

lS = 1 - 0.1*(log(f) - log(fc)).^2; 
lSc = max(lS);
S = exp(lS);


Sc = max(S);
reg = Sc/fc;

mindx = find(S > 0.5*Sc,1,'first');
maxdx = find(S > 0.5*Sc,1,'last');
fmin = f(mindx);
fmax = f(maxdx);

% no perturbation
Shat = S;
[Sc_,fc_,reg_,s,lS_] = fit_quadratic_density(f,Shat,fmin,fmax)

fdx = f>fmin & f<fmax;
fdx = find(fdx);
clf
loglog(f,exp(lS))
hold on
plot(f(fdx),exp(lS_))
vline(fmax)

% small perturbation
m = 1e3;
r = 0.1*lSc*randn(m,length(f));
Shat = exp(lS + r);
for idx=1:m
	[Sc_(idx,1),fc_(idx,1),reg_(idx,1),s_,S_] = fit_quadratic_density(f,Shat(idx,:),fmin,fmax);
	s.lfc(idx,1)     = s_.lfc;
	s.lSc(idx,1)     = s_.lSc;
	s.lreg(idx,1)    = s_.lreg;
	s.fc(idx,1)      = s_.fc;
	s.Sc(idx,1)      = s_.Sc;
	s.reg(idx,1)     = s_.reg;
	% raw is not a consistend estimator
	%[Sc_(idx,2),mdx] = max(Shat(idx,:));
	%fc_(idx,2)       = f(mdx);
	%reg_(idx,2)      = Sc_(idx,2)/fc_(idx,2);
end
plot(f,Shat(1,:))

'std(fc) fitted values, fit-est'
[std((fc_)),mean(s.fc)]

'std(Sc) fitted values, fit-est'
[std((Sc_)),mean(s.Sc)]

'std(reg) fitted values, fit est'
[std((reg_)),mean(s.reg)]

'E (fc),exact,fit'
[fc,mean((fc_))]

'E (Sc),exact,fit'
[Sc,mean((Sc_))]

'E (reg),exact,fit'
[reg,mean((reg_))]

if (0)
'E log(Sc)'
[lSc,mean(log(Sc_))]
'E reg, exact, fitted, raw'
[reg,mean(reg_)]
's(lreg)'
[std(log(reg_)),mean(s.lreg)]
%[std(log(Sc_)),mean(s.lSc)]
end
