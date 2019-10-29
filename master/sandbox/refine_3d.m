% Wed Jun 13 19:29:14 MSK 2012
% Karl Kästner, Berlin

% The Mathematics of Finite Elements and Applications X: MAFELAP 1999, p373

% consistency refinement type 1 - split a tetrahedron into 4-sub tetrahedra
function [P T Bc] = refine_3d(P, T, Bc, M)
	nt = size(T,1);
	np = size(P,1);
	% allocate memory
	T = [T; zeros(9*nt,4)];
	P = [P; zeros(nt,3)];
	Ph = java.util.Hashtable(6*size(T,1));

	for mdx=1:length(M)
	    tdx = M(mdx);
	    tetra_foursplit()
	
	% split domain boundary faces with 2D routine
	%TODO	
	end
	% resolve hanging nodes - 4 split
	%TODO	
	% resolve hanging nodes - 2 split
	%TODO	

	T = T(1:nt,:);
	P = P(1:np,:);
end % function refine_3d

% split a tetrahedron into four children
function [] = tetra_split_4(tdx,
	% remerge, if already unregularly refined
	% get the points
	p1 = T(tdx,1);
	p2 = T(tdx,2);
	p3 = T(tdx,3);
	p4 = T(tdx,4);

	% remove the old face references and mark neighbours for tetra_split_2
  	Fh = pop_face(Fh,p1,p2,p3,2);
	Fh = pop_face(Fh,p1,p2,p4,2);
	Fh = pop_face(Fh,p1,p3,p4,2);
	Fh = pop_face(Fh,p2,p3,p4,2);

	% find the face which has to be refined
	TODO permute such that p4 is opposit the side to be refined

	% four tetra-children
	T( tdx,:) = [p12 p13 p23 p4];
	T(np+1,:) = [ p1 p12 p13 p4];
	T(np+2,:) = [ p2 p12 p23 p4];
	T(np+3,:) = [ p3 p13 p23 p4];
	
	% faces of the four children
	Fh = push_faces(tdx,T,Fh);
	Fh = push_faces(np+1,T,Fh);
	Fh = push_faces(np+2,T,Fh);
	Fh = push_faces(np+3,T,Fh);
	np = np+3;
end % tetra_foursplit()

% split a tetrahedron into 8 children (regular refinement)
function [P np T nt Bc nb Ph Fh] = tetra_split_8(tdx,P,np,T,nt,Bc,nb,Ph,Fh)
	% remerge, if already unregularly refined once
	%   TODO

	% get the points
	    p1 = T(tdx,1);
	    p2 = T(tdx,2);
	    p3 = T(tdx,3);
	    p4 = T(tdx,4);

	% remove the old face references
  	   Fh = pop_face(Fh,p1,p2,p3);
	   Fh = pop_face(Fh,p1,p2,p4);
	   Fh = pop_face(Fh,p1,p3,p4);
	   Fh = pop_face(Fh,p2,p3,p4);

	% create six new side points
	    [P np Ph p12] = midpoint(P,np,Ph,p1,p2);
	    [P np Ph p13] = midpoint(P,np,Ph,p1,p3);
	    [P np Ph p14] = midpoint(P,np,Ph,p1,p4);
	    [P np Ph p23] = midpoint(P,np,Ph,p2,p3);
	    [P np Ph p24] = midpoint(P,np,Ph,p2,p4);
	    [P np Ph p34] = midpoint(P,np,Ph,p3,p4);

	% the four exterior tetrahedra
	   T( tdx,:) = [p1 p12 p13 p14];
           T(nt+1,:) = [p2 p12 p23 p24];
           T(nt+2,:) = [p3 p13 p23 p34];
           T(nt+3,:) = [p4 p14 p24 p34];
        % the four interior tetrahedra (three possebilities)
	% choose shortest edge as common edge
	d1234 = norm(P(p12,:)-P(p34,:));
	d1324 = norm(P(p13,:)-P(p24,:));
	d1423 = norm(P(p14,:)-P(p23,:));
	if (d1423 < d1324 && d1423 < d1234)
            T(nt+4,:) = [p12 p13 p14 p23];
	    T(nt+5,:) = [p12 p14 p23 p24];
	    T(nt+6,:) = [p13 p14 p23 p34];
            T(nt+7,:) = [p14 p23 p24 p34];
	elseif (d1234 < d1324)
	    T(nt+4,:) = [p12 p34 p13 p23];
	    T(nt+5,:) = [p12 p34 p13 p14];
	    T(nt+6,:) = [p12 p34 p14 p24];
	    T(nt+7,:) = [p12 p34 p23 p24];
	else
	    T(nt+4,:) = [p13 p24 p23 p34];
	    T(nt+5,:) = [p13 p24 p14 p34];
	    T(nt+6,:) = [p13 p24 p12 p23];
	    T(nt+7,:) = [p13 p24 p12 p14];
	end
	% push the new faces
	Fh = push_faces(tdx,T,Fh);
	for idx=nt:nt+7
		Fh = push_faces(nt+4,T,Fh);
	end
	nt = nt+7;
end % tetra_split_8


% splits a boundary face into four
function [Bc nb Ph] = split_boundary_4(idx, Bc, nb, Ph)
	p1 = Bc(idx,1);
	p2 = Bc(idx,2);
	p3 = Bc(idx,3);
	boundary = Bc(idx,4);
        p12 = Ph.get(hashkey(p1,p2));
        p13 = Ph.get(hashkey(p1,p3));
        p23 = Ph.get(hashkey(p2,p3));
	Bc(idx,:)  = [p12,p13,p23,boundary];
	Bc(nb+1,:) = [ p1,p12,p13,boundary];
	Bc(nb+2,:) = [ p2,p12,p23,boundary];
	Bc(nb+2,:) = [ p3,p13,p23,boundary];
	nb = nb+3;
end % split_boundary_4

% remove a face from the face->tetra relation hash
function [Fh Bc] = pop_face(Fh,Ph,Bc,nb,p1,p2,p3)
	key = haskey3(p1,p2,p3);
	val = Fh.remove(key);
	if (length(val) > 1)
		% other side not yet removed
		if (idx == val(1))
			val = val(2);
		else
			val = val(1);
		end
		if (val < 0)
			% boundary face
			[Bc nb Ph] = split_boundary_4(-val, Bc, nb, Ph);
		else
			TODO mark other tetra val for foursplit
			Fh.put(key,val);
		end
	end % else other side of face already split
end % pop_face

% push the four faces of a tetrahedron to the face->tetra relation hash
function Fh = push_faces(I,T,Fh)
	for idx=I
		key = hashkey3(T(idx,1),T(idx,2),T(idx,3));
		val = Fh.get(key);
		if (~isempty(val)) val(2) = idx; else val = idx; end
		Fh.put(hashkey3(key,val);

		key = hashkey3(T(idx,1),T(idx,2),T(idx,4));
		val = Fh.get(key);
		if (~isempty(val)) val(2) = idx; else val = idx; end
		Fh.put(hashkey3(key,val);

		key = hashkey3(T(idx,1),T(idx,3),T(idx,4));
		val = Fh.get(key);
		if (~isempty(val)) val(2) = idx; else val = idx; end
		Fh.put(hashkey3(key,val);

		key = hashkey3(T(idx,2),T(idx,3),T(idx,4));
		val = Fh.get(key);
		if (~isempty(val)) val(2) = idx; else val = idx; end
		Fh.put(hashkey3(key,val);
	end % for idx
end % push_faces()

function [P np Ph pdx] = midpoint(P,np,Ph,p1,p2)
 	% check side was already split and split if not
   	pdx = Ph.get(hashkey(p1,p2)); 
	if (isempty(pdx))
		np = np+1;
		P(np,1) = 0.5*(P(p1,1) + P(p2,1));
		P(np,2) = 0.5*(P(p1,2) + P(p2,2));
		P(np,3) = 0.5*(P(p1,3) + P(p2,3));
		Ph.put(hashkey(p1,p2),np);
		pdx = np;
	end
end % midpoint()

%{
	%
	%
	%

% splits a line into two segments
% P : array of point coordinates
% L : array of lines, entries are indices into P
function [P np L nl] = split_1d_2(ldx, P, np, L, nl)
	% get old point indices
	p1 = L(ldx,1);
	p2 = L(ldx,2);
	% check that the side has not yet been split
	todo
	% new interior point
	np = np+1;
	P(np,1) = 0.5*(P(p1,1) + P(p2,1));
	P(np,2) = 0.5*(P(p1,2) + P(p2,2));
	P(np,3) = 0.5*(P(p1,3) + P(p2,3));
	% update old line segment
	L(ldx,2) = np;
	% new line segment
	nl = nl+1;
	L(nl,1) = np;
	L(nl,2) = p2;
end % split_1d_2()

% splits a triangular face into four triangles
function [P np F nf L nl] = split_2d_4(fdx,P,np,L,nl,F,nf)
	% check, that this triangle has not yet been split
	todo, 
	% split the three sides
	for idx=1:3
		[P np L nl] = split_1d_2(F(fdx,idx), P, np, L, nl);
	end
	% update interior triangle
	% three new corner children
end

% regular refinement - splits a tetrahedron into 8 tetrahedra
% tdx : index of tetrahedron to be splitted
function [P, np, L, nl, F, nf, T, nt] = split_3d_8(tdx, P, np, L, nl, F, nf, T, nt)
	% split the four faces
	for idx=1:4
		[P np L nl F nf] = split_2d_4(T(tdx,idx),P,np,L,nl,F,nf);
	end
	% four corner children
	% four interior children
end

%}

