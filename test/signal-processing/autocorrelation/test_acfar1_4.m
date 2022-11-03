% 2019-06-27 18:24:10.402420751 +0200
r = [1-2.^-(0:20)];  clf; for idx=1:length(r); plot(acfar1(r(idx),100,0:100)); hold on; end; legend(num2str(cvec(r)))

