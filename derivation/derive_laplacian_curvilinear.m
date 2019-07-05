% Tue 13 Mar 13:06:38 CET 2018
% Wed 29 Aug 09:24:06 CEST 2018
% c.f Tryggvason 2011, Computational Fluid Dynamics Boundary-Fitted Coordinates
% note: trygvasson gives also conservative forms

reset(symengine)

if (1)

% syms f(x) x(chi) chi, f = f(x), diff(f,chi), functionalDerivative(f,x) 
syms dx_ds ds_dx df_ds

disp('1d')
ds_dx = solve(dx_ds*ds_dx == 1,ds_dx)
df_dx = df_ds*ds_dx

disp('2d')
syms df_dn dx_dn dn_dx dy_ds dy_dn ds_dx df_dx df_dy
s = solve([dx_ds*ds_dx + dx_dn*dn_dx == 1, 
           dy_ds*ds_dx + dy_dn*dn_dx == 1],[ds_dx,dn_dx]);
ds_dx  = s.ds_dx
dn_dx  = s.dn_dx
%df_dx = df_ds*ds_dx + df_dn*dn_dx
%s=solve([df_dx == df_ds*ds_dx + df_dn*dn_dx
s=solve([df_ds == df_dx*dx_ds + df_dy*dy_ds;
         df_dn == df_dx*dx_dn + df_dy*dy_dn],[df_dx,df_dy]);
df_dx = s.df_dx
df_dy = s.df_dy

% ---

syms f(s,n) x(s,n) y(s,n) J_
%J(s,n)

J = diff(x,s)*diff(y,n) - diff(x,n)*diff(y,s)

df_dx = @(f) 1/J*(diff(y,n)*diff(f,s) - diff(y,s)*diff(f,n))
df_dy = @(f) 1/J*(diff(x,s)*diff(f,n) - diff(x,n)*diff(f,s))


d2f_dx2 = @(f) df_dx(df_dx(f))
d2f_dxy = @(f) df_dx(df_dy(f))
d2f_dy2 = @(f) df_dy(df_dy(f))

%df2_dx2_ = d2f_dx2(f)
%d2f_dx2_ = d2f_dx2(f)

d = {df_dx(f), df_dy(f), d2f_dx2(f), d2f_dxy(f), d2f_dy2(f)}

% make human readable
for idx=1:length(d)

syms y_ss y_sn y_nn y_s y_n
syms x_ss x_sn x_nn x_s x_n
syms f_ss f_sn f_nn f_s f_n

d{idx} = simplify(d{idx})

% first derivatives
d{idx} = subs(d{idx}, diff(y,s), y_s);
d{idx} = subs(d{idx}, diff(y,n), y_n);
d{idx} = subs(d{idx}, diff(x,s), x_s);
d{idx} = subs(d{idx}, diff(x,n), x_n);
d{idx} = subs(d{idx}, diff(f,s), f_s);
d{idx} = subs(d{idx}, diff(f,n), f_n);
% second derivatives
d{idx} = subs(d{idx}, diff(x,s,2), x_ss);
d{idx} = subs(d{idx}, diff(x,n,2), x_nn);
d{idx} = subs(d{idx}, diff(y,s,2), y_ss);
d{idx} = subs(d{idx}, diff(y,n,2), y_nn);
d{idx} = subs(d{idx}, diff(diff(y,s),n), y_sn);
d{idx} = subs(d{idx}, diff(diff(x,s),n), x_sn);
d{idx} = subs(d{idx}, diff(f,s,2), f_ss);
d{idx} = subs(d{idx}, diff(f,n,2), f_nn);
d{idx} = subs(d{idx}, diff(diff(f,s),n), f_sn);

%d{idx} = subs(d{idx},x_n*y_s - x_s*y_n,J_)

end % for idx
df_dx = d{1}
df_dy = d{2}
d2f_dx2 = d{3}
d2f_dxy = d{4}
d2f_dy2 = d{5}
L = simplify(d2f_dx2+d2f_dy2)

d2f_dx2 = collect(collect(collect(collect(collect(d2f_dx2,f_ss),f_sn),f_nn),f_s),f_n)
d2f_dxy = collect(collect(collect(collect(collect(d2f_dxy,f_ss),f_sn),f_nn),f_s),f_n)
d2f_dy2 = collect(collect(collect(collect(collect(d2f_dy2,f_ss),f_sn),f_nn),f_s),f_n)
L       = collect(collect(collect(collect(collect(L,f_ss),f_sn),f_nn),f_s),f_n)

%d2f_dx2_ = d{1}
%d2f_dy2_ = d{2}
%L    

if (1)
	df_dx = (subs(df_dx,x_s*y_n - x_n*y_s,J_))
	df_dy = (subs(df_dy,x_s*y_n - x_n*y_s,J_))
	d2f_dx2 = (subs(d2f_dx2,x_s*y_n - x_n*y_s,J_))
	d2f_dxy = (subs(d2f_dxy,x_s*y_n - x_n*y_s,J_))
	d2f_dy2 = (subs(d2f_dy2,x_s*y_n - x_n*y_s,J_))
	L       = (subs(L,x_s*y_n - x_n*y_s,J_))
end

end

% with vandermonde
if (0)
row = ['l','c','r']
col = ['t','m','b']
syms f0 df_dx df_dy d2f_dx2 d2f_dxy d2f_dy2
f_ = [f0 df_dx df_dy d2f_dx2 d2f_dxy d2f_dy2].'
syms xcm ycm A A_ b  AA
for idx=1:3
 for jdx=1:3
	eval(['syms x',row(idx),col(jdx)])
	eval(['syms y',row(idx),col(jdx)])
	eval(['syms f',row(idx),col(jdx)])
	x = eval(['x',row(idx),col(jdx)])
	y = eval(['y',row(idx),col(jdx)])
	f = eval(['f',row(idx),col(jdx)])
	A_(idx+(jdx-1)*3,1) = f == vander_2d((x-xcm),(y-ycm),2)*f_
	eval(['syms dx',row(idx),col(jdx)])
	eval(['syms dy',row(idx),col(jdx)])
	eval(['syms f',row(idx),col(jdx)])
	dx = eval(['dx',row(idx),col(jdx)])
	dy = eval(['dy',row(idx),col(jdx)])
	if (idx==2 && jdx == 2)
		dx=0; dy=0;
	end
	f = eval(['f',row(idx),col(jdx)])
	A(idx+(jdx-1)*3,1:6) = vander_2d(dx,dy,2)
	b(idx+(jdx-1)*3,1) = f
 end
end
%df = simplify(A.'*A) \ simplify(A.'*b)
for idx=1:9
for jdx=1:6
	eval(['syms a',num2str(idx),num2str(jdx)])
	a = eval(['a',num2str(idx),num2str(jdx)])
	AA(idx,jdx) = a
end
end
tic(); iA=inv(AA.'*AA)*AA.',toc
end
