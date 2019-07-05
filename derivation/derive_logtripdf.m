% Fri 23 Mar 12:24:28 CET 2018

syms x a b c den F x positive

fl = @(x) (log(x/a)/x)/(log(b/a)/b)*1/den
fr = @(x) (log(c/x)/x)/(log(c/b)/b)*1/den
if (0)
	% alternative (linear pdf, but not linear histogram
	fl = @(x)(log(x) - log(a))/(log(b)-log(a))*1/den
	fr = @(x) (log(c) - log(x))/(log(c)-log(b))*1/den
end

Fl = int(fl,x)
Fl = Fl - subs(Fl,x,a)
Fl = simplify(Fl,'ignoreanalyticconstraints',true)

Fr = int(fr,x)
Fr = Fr - subs(Fr,x,b) + subs(Fl,x,b)
Fr = simplify(Fr,'ignoreanalyticconstraints',true)

%den = subs(Fl,x,b) + subs(Fr,x,c)
den_ = subs(Fr,x,c)
den_ = simplify(den_);

%fl = @(x) fl(x)/den
%fr = @(x) fr(x)/den
%end


Fb = subs(Fl,b)
%subs(Fl,x,b) - subs(Fl,x,a);
Fb = simplify(Fb)


iFl = solve(Fl-F,x)
iFr = solve(Fr-F,x)

%iFl = solve(Fl-F,x)
%iFl = -F/wrightOmega(-log(-e/(F*a*den*log(a/b))))*(-den*(log(a) - log(b)))

