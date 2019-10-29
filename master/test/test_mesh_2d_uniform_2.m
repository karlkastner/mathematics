javaaddpath('.')
javaaddpath('/usr/share/java/jama.jar');
figure(1); clf

subplot(2,2,1)
n = [5 5], L = [1 1]
[P T Bc] = mesh_2d_uniform(n, L);
display_2d(P,T,Bc,8);
axis equal
axis tight

subplot(2,2,2)
n = [10 5], L = [2 1]
[P T Bc] = mesh_2d_uniform(n, L);
display_2d(P,T,Bc,8);
axis equal
axis tight

subplot(2,2,3)
n = [5 10], L = [1 2]
[P T Bc] = mesh_2d_uniform(n, L);
display_2d(P,T,Bc,8);
axis equal
axis tight

subplot(2,2,4)
X = (0:5)'/5;
Y = (0:10)'/5;
[P T Bc] = mesh_2d_uniform([], {X,Y});
display_2d(P,T,Bc,8+2);
axis equal
axis tight

B = java_assemble_2d_phi_phi(P, T, [], @int_2d_gauss_3);
figure(2); clf
spy(B);

