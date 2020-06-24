function [ p ] = multi_randperm( n,k )

if k<factorial(n)
    p(1,:)=1:n;
    for i=2:k
        p(i,:)=randperm(n);
    end
else
    p=perms(1:n);
    p=p(end:-1:1,:);
end
end

