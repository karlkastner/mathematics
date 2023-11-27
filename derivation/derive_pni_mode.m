% 2023-10-27 10:00:40.057038956 +0200
% Karl Kästner, Berlin
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
syms f a1 a2;
		 S = pni_pdf(f,a1,a2);
		 dS_df = diff(S,f);
		 [nS,dS] = numden(dS_df);
		 nS=simplify(nS);
		 f0 = solve(nS==0,f)
		Sc = subs(S,f,f0(3))

 syms a1 a2 real;
		 iS=simplify(iS,'ignoreanalyticconstraints',true);
		 liS = limit(iS,f,inf)

		[fc,Sc]=pni_mode(a1,a2)
syms fc0 Sc0
		a10 = solve(fc-fc0,a1);
		eq = subs(Sc-Sc0,a1,a10);
		[n,d] = numden(eq)


%(2^(1/2)*signIm((2*a2*(4*a2 - a1^2)^(3/2) + a1*a2^2*16i - a1^3*a2*8i - a1^2*(4*a2 - a1^2)^(3/2) + a1^5*1i)/((4*a2 - a1^2)*(8*a2^2 - 6*a1^2*a2 + a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2)))*(8*a2^2 - 6*a1^2*a2 + a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2))/(8*a1*(4*a2 - a1^2)) - (2^(1/2)*signIm(((a1^2*(4*a2 - a1^2)^(3/2) + a1*a2^2*16i - a1^3*a2*8i - 2*a2*(4*a2 - a1^2)^(3/2) + a1^5*1i)*1i)/((4*a2 - a1^2)*(6*a1^2*a2 - 8*a2^2 - a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2)))*(6*a1^2*a2 - 8*a2^2 - a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2)*1i)/(8*a1*(4*a2 - a1^2))
%signIm((2*a2*(4*a2 - a1^2)^(3/2) + a1*a2^2*16i - a1^3*a2*8i - a1^2*(4*a2 - a1^2)^(3/2) + a1^5*1i)/((4*a2 - a1^2)*(8*a2^2 - 6*a1^2*a2 + a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2)))
%signIm(((a1^2*(4*a2 - a1^2)^(3/2) + a1*a2^2*16i - a1^3*a2*8i - 2*a2*(4*a2 - a1^2)^(3/2) + a1^5*1i)*1i)/((4*a2 - a1^2)*(6*a1^2*a2 - 8*a2^2 - a1^4 + a1*(4*a2 - a1^2)^(3/2)*1i)^(1/2)))

pi = sym(pi);

eq = Sc0 - 4*2^(1/2)*(- 4*a2^2*fc0^2*pi^2 + a2)^(1/2) - 16*Sc0*a2^2*fc0^4*pi^4
eq = (Sc0 - 16*Sc0*a2^2*fc0^4*pi^4)^2 - (4*2^(1/2)*(- 4*a2^2*fc0^2*pi^2 + a2)^(1/2))^2
