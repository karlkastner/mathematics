% 2023-10-27 10:01:22.001736235 +0200
		a1 = 0.01
		a2 = 0.05
		[fc,Sc]=pni_mode(a1,a2);
          
		%f0 = matlabFunction(f0(3));
		%dS = matlabFunction(dS_df);
		 f = linspace(0,10,1e3)';
		 S=pni_pdf(f,a1,a2);
		 plot(f,S);
		%S,cdiff(S)./cdiff(f),dS(a1,a2,f)]);
		vline(fc);
		hline(Sc)
		 hline(0);
%		 vline(f0(a1,a2))


		[a1,a2] = pni_mode2par(fc,Sc)
