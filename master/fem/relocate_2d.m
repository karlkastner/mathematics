% Thu Jun 14 12:12:06 MSK 2012
% Karl Kästner, Berlin

%	r-convergence: cauchy's method with line search
 %       - highly nonlinear optimisation problem
  %      - gradient and line search calculation can be reduced by assuming that a point loction
   %             mainly effects the contributions by neighbouring triangles
    %    - todo: recalculate dV
    %            for each triangle
    %                    check if contained in the triangle
    %                            interpolate value
    %                            break
    %    - derivative calculation
    %            - w = 0.5, if either dx or dy is zero, set to 1
    %            - only weight is w for all entries


function P_ = relocate_2d(P, T, Bc, N, err, nH, dV, area, h_side, s_angle, C)
	lambda = 1;
	lp = size(P,1);
	df = zeros(size(P)); 

	% get an inverse point->triangle map
	Np = zeros(size(P,1),7);
	Np(:,1) = 1;
	for tdx=1:size(T,1)
		for idx=1:3 % todo, four in 3D
			Np(T(tdx,idx),1) = Np(T(tdx,idx),1) + 1;
			Np(T(tdx,idx),Np(T(tdx,idx),1)) = tdx;
		end
	end

	% calculate the gradient
	for pdx=1:lp
		% choose dx,dy according to minimum side length of neighbouring triangles
		h_min = min(min(h_side(Np(pdx,2:Np(pdx,1)))));
		dx = 0.1*h_min;
		% save point value
		P_backup = P(pdx,:);
		err_old = sum(err(Np(pdx,2:Np(pdx,1))));
		% for each dimension
		for idx=1:size(P,2)
			% change point value
			P(pdx,idx) = P(pdx,idx) + dx;
			% TODO interpolate the value at the point and recalculate dV
			% recalculate secondary triangle properties of neighbouring elements
			[area h_side s_angle C] = recalculate_regularity_2d(Np(pdx,2:Np(pdx,1)), P,T,area,h_side,s_angle,C);
			% recalculate the error contribution of neighbouring elements
			[err nH] = estimate_err_2d_3(Np(pdx,2:Np(pdx,1)), N, dV, area, h_side, s_angle, C, err, nH);
			%err_new = sum(err(Np(pdx,2:Np(pdx,1))));
			err_p = sum(err(Np(pdx,2:Np(pdx,1))));
			% restore the point value
			P(pdx,:) = P_backup;

			% change point value
			P(pdx,idx) = P(pdx,idx) - dx;
			% TODO interpolate the value at the point and recalculate dV
			% recalculate secondary triangle properties of neighbouring elements
			[area h_side s_angle C] = recalculate_regularity_2d(Np(pdx,2:Np(pdx,1)), P,T,area,h_side,s_angle,C);
			% recalculate the error contribution of neighbouring elements
			[err nH] = estimate_err_2d_3(Np(pdx,2:Np(pdx,1)), N, dV, area, h_side, s_angle, C, err, nH);
			%err_new = sum(err(Np(pdx,2:Np(pdx,1))));
			err_m = sum(err(Np(pdx,2:Np(pdx,1))));
			% restore the point value
			P(pdx,:) = P_backup;

			% get the gradient
			df(pdx,idx) = -(err_p - err_m)/(2*dx);
			%df(pdx,idx) = (err_new - err_old)/dx;
		end % for idx
		% line search into gradient direction
		ndf = norm(df(pdx,:));
		if (ndf > 1e-12*h_min)
			d = 0.5*h_min*df(pdx,:)/ndf;
			for jdx=1:5
				% change point value
				P(pdx,:) = P_backup + d;
				% recalculate secondary triangle properties of neighbouring elements
				[area h_side s_angle C] = recalculate_regularity_2d(Np(pdx,2:Np(pdx,1)), P,T,area,h_side,s_angle,C);
				% recalculate the error contribution of neighbouring elements
				[err nH] = estimate_err_2d_3(Np(pdx,2:Np(pdx,1)), N, dV, area, h_side, s_angle, C, err, nH);
				%err_new = sum(err(Np(pdx,2:Np(pdx,1))));
				err_new = sum(err(Np(pdx,2:Np(pdx,1))));
				if (err_new < err_old) break; end
				% reduce step width
				d = 0.5*d;
			end
			df(pdx,:) = d;
		end
		% restore the point value
		P(pdx,:) = P_backup;
	end % for pdx

	% relocate the points by one cauchy step
	P_ = lambda*df + P;

	% restore locations of boundary points
	% TODO: use a back-projection (take care of corners)
	for idx=1:size(Bc,1)
		P_(Bc(idx,1),:) = P(Bc(idx,1),:);
		P_(Bc(idx,2),:) = P(Bc(idx,2),:);
	end
	% restore the point at zero
	[min_ fdx] = min(sum(abs(P),2));
	P_(fdx,:) = P(fdx,:);
end % function relocate_2d

%{
function P_ = relocate_2d(P,T,Bc,nH,h_side,s_angle,C)
	% allocate memory
	P_ = zeros(size(P));
	w = zeros(size(P));
	% contributions to new point coordinates (todo vectorise)
	for tdx=1:size(T,1)
		for idx=1:3
			%w_ti = nH(tdx) * (max(h_side(tdx,:))/s_angle(tdx,idx))^2;
			w_ti = nH(tdx) * (h_side(tdx,idx)/s_angle(tdx,idx))^2;
			%w_ti = err(tdx);
			P_(T(tdx,idx),:) = P_(T(tdx,idx),:) + C(tdx,:) * w_ti;
			w(T(tdx,idx),:) =  w(T(tdx,idx),:) + w_ti;
		end % for idx
	end % for tdx
	% norm new point coordinates
	P_ = P_./w;
	% restore locations of boundary points
	% TODO: use a back-projection (take care of corners)
	for idx=1:size(Bc,1)
		P_(Bc(idx,1),:) = P(Bc(idx,1),:);
		P_(Bc(idx,2),:) = P(Bc(idx,2),:);
	end
end % function relocate_2d
%}

