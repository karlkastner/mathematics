% Fri Jul  4 14:45:50 WIB 2014
% Karl Kastner, Berlin
%
%% use java port for speed up
%
function [S N pid sid] = xy2sn_java(cX,cY,cS,cSeg,X,Y)
	% there is a class name conflict
	%javarmpath('/home/pia/Documents/master/thesis/src/fem/class');
	try
		m = size(X,2);
		if (m > 1)
			X = flat(X);
			Y = flat(Y);
			%S = flat(S);
			%N = flat(N);
		end

		XY2sn_java = javaObject('XY2sn');
		cX = double(cX);
		cY = double(cY);
		X  = double(X);
		Y  = double(Y);
		cSeg = int32(cSeg);
		S  = zeros(size(X));
		N  = zeros(size(X));
		XY2sn_java.xy2sn_(cX,cY,cS,cSeg,X,Y,S,N);

	catch e
		e
		e.message
		error('Error in xy2sn_java');
	end
	S   = XY2sn_java.S;
	N   = XY2sn_java.N;
	pid  = double(XY2sn_java.Nkey);
	if (m > 1)
		S   = reshape(S,[],m);	
		N   = reshape(N,[],m);	
		pid  = reshape(pid,[],m);
	end
	% quick fix for bug for first element in qtree.nn 
	pid(pid == 0) = 1;
	sid = cSeg(pid);
	% restore the path
%	javaaddath('/home/pia/Documents/master/thesis/src/fem/class');
end % xy2sn_java

