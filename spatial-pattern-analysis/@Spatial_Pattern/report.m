% Wed  8 Feb 12:17:32 CET 2023
function report(obj)
	lambda_c = obj.lambda_c();
	if (obj.stat.isisotropic)
		fprintf('Isotropic\n');
		fprintf('S_{cr}/lambda_c  %f\n',obj.stat.Sc.radial.clip/lambda_c);
		fprintf('S_{ct}           %f\n',obj.stat.Sc.angular.clip);
		fprintf('L_{eff}/lambda_c %f\n',obj.stat.L_eff.r/lambda_c);
	else
		fprintf('Anisotropic\n');
		fprintf('S_{cx}/lambda_c  %f\n',obj.stat.Sc.x.clip/lambda_c);
		fprintf('S_{cy}/lambda_c  %f\n',obj.stat.Sc.y.clip/lambda_c);
		fprintf('L_{eff}/lambda_c %f\n',obj.stat.L_eff.r/lambda_c);
	end
	fprintf('p_{periodic} %f\n',obj.stat.p_periodic);
end

