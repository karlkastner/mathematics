clear
aa = rand(3,3,1);


eig_ = eig3(aa);

for idx=1:size(aa,3)
	eig__(idx,:) = eig(aa(:,:,idx));
end

eig_ = sort(eig_')';
eig__ = sort(eig__')';

eig_
eig__
eig_-eig__
