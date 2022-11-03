
% n.b. matlab normalized the incomplete beta function
beta(a,b)*betainc(x,a,b)
mybetainc(x,a,b,n)
% identity
(-1)^a*mybetainc(x./(x-1),a,1-a-b,n)
