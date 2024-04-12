% Thu  8 Feb 12:31:14 CET 2024
function znext = step_bdf2(I,A,z,dt)
	znext = (I - (2*dt/3)*A) \ ((4/3)*z(:,1) - (1/3)*z(:,2));
end

