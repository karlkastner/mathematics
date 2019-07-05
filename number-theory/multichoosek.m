% 2018-10-11 13:25:27.874645589 +0200
% Karl Kastner, Berlin
% implementation by knedlsepp Feb 2 2015 on stackoverflow 
%
%% all combinations of lenght k from set values with repetitions
%% c.f. nchoosek, combinations without repetition
%%
%% input :
%% 	x : scalar integer or vector of arbitrary numbers
%%	k : length of subsets
%% output :
%%	if x scalar : number of combinations  
%%	if x vector : the exact combinations
%%
function combs = nmultichoosek(values, k)
if numel(values)==1 
    n = values;
    combs = nchoosek(n+k-1,k);
else
    n     = numel(values);
    combs = bsxfun(@minus, nchoosek(1:n+k-1,k), 0:k-1);
    combs = reshape(values(combs),[],k);
end
