% Wed Jul 11 17:25:17 MSK 2012
% Karl Kästner, Berlin

% coordinates and weights for numerical Gauss quadrature
% 7th-order accurate
function [w, b, flag] = int_1d_gauss_4()
	 w = [0.326072577431273;
	      0.326072577431273;
              0.173927422568727;
              0.173927422568727];
	 b = [0.330009478207572 0.669990521792428;
              0.669990521792428 0.330009478207572;
	      0.069431844202974 0.930568155797026;
              0.930568155797026 0.069431844202974];
	flag = 0;
end % int_2d_gauss_4

