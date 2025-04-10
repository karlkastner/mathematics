% Wed  8 Feb 12:17:32 CET 2023
% Karl Kästner, Berlin
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%% report statistics of analysis
%
function report(obj)
	lambda_c = obj.lambda_c();
	field = obj.opt.scalefield;
	if (obj.stat.isisotropic)
		fprintf('Isotropic\n');
		fprintf('S_{rc}/lambda_c  %f\n',obj.stat.Sc.radial.(field)/lambda_c);
		fprintf('S_{tc}^+         %f\n',obj.stat.Sc.angular_p.(field));
		fprintf('L_{eff}/lambda_c %f\n',obj.stat.L_eff.r/lambda_c);
	else
		fprintf('Anisotropic\n');
		fprintf('S_{xc}^+/lambda_c  %f\n',obj.stat.Sc.xp.(field)/lambda_c);
		fprintf('S_{yc}/lambda_c  %f\n',obj.stat.Sc.y.(field)/lambda_c);
		fprintf('L_{eff}/lambda_c %f\n',obj.stat.L_eff.r/lambda_c);
	end
	fprintf('\lambda_c        %f\n',lambda_c);
	fprintf('p_{periodic} %f\n',obj.stat.p_periodic);
	fprintf('L %f %f\n',obj.L);
end

