x = rand(1e6,3); y =zeros(1e6,2); tic; y(:,1)=gampdf(x(:,1),x(:,2),x(:,3)); toc, tic; y(:,2)=gampdf_man(x(:,1),x(:,2),x(:,3)); toc, rms(y(:,2)-y(:,1))
