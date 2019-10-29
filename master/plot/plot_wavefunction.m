% Sun Nov 27 22:20:02 MSK 2011
% Karl Kästner, Berlin

function plot_wavefunction(X,V,dimension)
 	%E = fliplr(E(:,end)');
	 %plot(X, W*N(1)*norm(E) + W.^0*diag(E),'k');
	n = length(X);
	switch (dimension)
		case {1}
			plot(X,V);
 		case {2}
			V2 = reshape(V,n,n);
			imagesc(V2);
			axis square
			axis off
 		case {3}
			V3 = reshape(V,n,n,n);
			lim = 0.5*max(V);
			[XX YY ZZ] = meshgrid(X,X,X);
			p = patch(isosurface(XX,YY,ZZ,V3,lim));
			%isonormals(XX,YY,ZZ,W3_,p)
			set(p,'FaceColor','red','EdgeColor','none');
			daspect([1,1,1])
			view(3);
			camlight 
			lighting gouraud
			axis([X(1) X(end) X(1) X(end) X(1) X(end)])
			grid on
 	end % switch
end % function plot_wavefunction

