function [ d,d2 ] = permutation_dimension( perm )
%PERMUTATION_DIMENSION determines the Haussdorff dimension of the
%permutation fractal

eps=1/(length(perm));
len=sum(abs(diff(perm)));
d=-log(len)/log(eps);
len2=sum(sqrt(diff(perm).^2+1));
d2=-log(len2)/log(eps);

%Formula:
%s=log2(length(perm));
%dim=@(x,m) -log(sum(2.^(m-x).*sqrt(2.^(2*x-2)+1)))/log(2^(-m));
%dim(1:s,s)

%Or even better
%dim_adv=@(x,m) 1+log(sum(2.^-x.*sqrt(2.^(2*x-2)+1)))/(m*log(2));
end