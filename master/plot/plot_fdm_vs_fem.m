f = {'../dat-new/fdm-1338379748.mat', '../dat-new/fdm-1338382230.mat', ...
	'../dat-new/fem-2d-1338937673.mat', '../dat-new/fem-2d-1338938032.mat'}
%	'../dat-new/fem-2d-1338300608.mat', '../dat-new/fem-2d-1338301105.mat'}
colour = {'b', 'g', 'r', 'm'}

pattern = {'.--b', '.-b', '.--r', '.-r' }

figure(1); clf; subplot(1.5,1.5,1)


for idx=1:4
	load(f{idx});
	
	        err_true = abs ( ( E - E_true*ones(1,size(E,2))).^2 ...
				./ (E_true*ones(1,size(E,2))));
		nerr_true = sqrt(sum(err_true,1));
	%Err = (E - E_true*ones(1,length(N(:,2))))./(E_true*ones(1,length(N(:,2))));
	%nErr = sqrt(sum(Err.^2,1));

	%loglog(N(:,1), err_est,'.-','color',colour{idx});
	if (idx < 3)
		%loglog(cumsum(sum(T)), nerr_true,'.-','color',colour{idx});
		loglog(cumsum(sum(T)), nerr_true,pattern{idx});
	else
		%loglog(cumsum(sum(Ti,2)), nerr_true,'.-','color',colour{idx});
		loglog(cumsum(sum(Tr,2)), nerr_true,pattern{idx});
	end
	grid on; set(gca,'minorgrid','none');
	title('Effciency');
	xlabel('time [s] (accumulated)')
	ylabel('true |\lambda - \lambda_*|/|\lambda_*|');
	hold on;
	xlim([1e-1 1e3])
	ylim([5e-5 5e0])
end
legend('fdm k = 1', 'fdm k = 23', 'fem k =  1', 'fem k = 23');

print('-depsc', '../dat-new/fdm-vs-fem-2d-efficency');

