% Fri  4 Oct 14:57:13 CEST 2024
	syms z zl  bl br zr real
	syms C h hl hr dx positive

	Sl = (h + z  - hl - zl)/dx;
	Sr = (hr+ zr - h  - z )/dx;

	% velocity, this is the coefficient
	ul = -C*sign(Sl)*sqrt(0.5*(hl + h ).*abs(Sl)) 
	ur = -C*sign(Sr)*sqrt(0.5*(h  + hr).*abs(Sr)) 

Jl = diff(dh_dt,hl)
Jr = diff(dh_dt,hr)
Jc = diff(dh_dt,h )

	dul_dhl = diff(ul,hl)
	dur_dhr = diff(ur,hr)
	dul_dhc = diff(ul,h)
	dur_dhc = diff(ur,h)


	hli = 0.5*(hl + h);
	hri = 0.5*(hr + h);
	dul_dhl_ =  ul.*(0.25./hli - 0.5./(dx*Sl))
	dur_dhr_ =  ur.*(0.25./hri + 0.5./(dx*Sr))

	dul_dhc_ =  ul.*(0.25./hli + 0.5./(dx*Sl))
	dur_dhc_ =  ur.*(0.25./hri - 0.5./(dx*Sr))

	dh_dt = -(    0.5*ur.*( (1-br).*h   + (1+br).*hr) ...
	            - 0.5*ul.*( (1-bl).*hl  + (1+bl).*h ) ...
		 )/dx;


	%$dur_dhr = diff(ur,hr);
	Jl_  = (   0.5*ul.*(1-bl) ...
		      +0.5*dul_dhl_.*( (1-bl).*hl  + (1+bl).*h ) ...
		    )/dx;
	Jr_  = -(   0.5*ur.*(1+br) ...
		   +0.5*dur_dhr_*( (1-br).*h   + (1+br).*hr) ...
	        )/dx;

	Jc_ =  -( ...
		    0.5*ur.*(1-br) ...
		   +0.5*dur_dhc_.*((1-br).*h   + (1+br).*hr) ...
		   -0.5*ul.*(1+bl) ...
		   -0.5*dul_dhc_.*((1-bl).*hl  + (1+bl).*h ) ...
	       )/dx;
	Jc__ = -(Jl_ + Jr_);

du_dhl    = diff(ul,hl)
f.dul_dhl  = matlabFunction(dul_dhl);
f_.dul_dhl = matlabFunction(dul_dhl_);
f.dur_dhr  = matlabFunction(dur_dhr);
f_.dur_dhr = matlabFunction(dur_dhr_);
f.dur_dhc  = matlabFunction(dur_dhc);
f_.dur_dhc = matlabFunction(dur_dhc_);
f.Jl      = matlabFunction(Jl);
f_.Jl      = matlabFunction(Jl_);
f.Jr      = matlabFunction(Jr);
f_.Jr      = matlabFunction(Jr_);
f.Jc      = matlabFunction(Jc);
f_.Jc      = matlabFunction(Jc_);
f__.Jc      = matlabFunction(Jc__);
if (0)
f.dul_dhl(1,1,1,2,0,0)
f_.dul_dhl(1,1,1,2,0,0)  	
f.dul_dhl(1,1,2,1,0,0) 
f_.dul_dhl(1,1,2,1,0,0) 
'dur_dhr'
f.dur_dhr(1,1,1,2,0,0)
f_.dur_dhr(1,1,1,2,0,0)  	
f.dur_dhr(1,1,2,1,0,0) 
f_.dur_dhr(1,1,2,1,0,0) 
end
f.dur_dhc(1,1,1,2,0,0)
f_.dur_dhc(1,1,1,2,0,0)  	
f.dur_dhc(1,1,2,1,0,0) 
f_.dur_dhc(1,1,2,1,0,0) 
pause
%    (C,bl,br,dx,h,hl,hr,z,zl,zr)
f.Jc( 1, 0, 0, 1,1, 2, 3,0,0,0)
f_.Jc(1, 0, 0, 1,1, 2, 3,0,0,0)
f__.Jc(1, 0, 0, 1,1, 2, 3,0,0,0)
%    (C,bl,br,dx,h,hl,hr,z,zl,zr)
f.Jc( 1, 0, 0, 1,1, 2, 0.5,0,0,0)
f_.Jc(1, 0, 0, 1,1, 2, 0.5,0,0,0)
f__.Jc(1, 0, 0, 1,1, 2, 0.5,0,0,0)

f.Jc( 1, 0, 0, 1,2, 1, 0.5,0,0,0)
f_.Jc(1, 0, 0, 1,2, 1, 0.5,0,0,0)
f__.Jc(1, 0, 0, 1,2, 1, 0.5,0,0,0)

