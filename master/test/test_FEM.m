javaaddpath('.')
javaaddpath('/usr/share/java/jama.jar');

n = 2; L0 = [1 1];

int = @int_2d_gauss_3;
[P T Bc] = mesh_2d_uniform(n, L0);
%int = @int_2d_gauss_3;
%[P T Bc] = promote_2d_3_6(P,T,Bc);
int = @int_2d_gauss_6;
[P T Bc] = promote_2d_3_10(P,T,Bc);
%[P T Bc] = promote_2d_3_15(P,T,Bc);
%int = @int_2d_gauss_12;


'testing phi phi'
tic();
[A_ buf_] = assemble_2d_phi_phi(P, T, @f_coulomb, int);
toc()

tic();
[w b flag] = feval(int); 
[A buf] = java_assemble_2d_phi_phi(P, T, Func_Coulomb, int);
toc()
%pause
%[buf_(:,end) buf(:,end) buf_(:,end)-buf(:,end)], full(A), full(A_)
sum(sum((A_ - A).^2))

'testing dphi dphi'
tic();
[A_ buf_] = assemble_2d_dphi_dphi(P, T, [], int);
toc()
tic();
[A buf] = java_assemble_2d_dphi_dphi(P, T, [], int);
toc()
%[buf_(:,end) buf(:,end) buf_(:,end)-buf(:,end)]
% full(A), full(A_)
sum(sum((A_ - A).^2))

