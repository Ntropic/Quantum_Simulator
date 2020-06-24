function [ x,n2 ] = Order_Of_Convergence( D,n )
%Calculates the order of convergence
eps=10^-15;
if nargin==2
    n2=n(2:end);
    x=zeros(1,length(n2));
    
    for i=1:length(n2)
        if D(i+1)>eps && D(i)>eps
            x(i)=(log(D(i+1))-log(D(i)))/log((n(i))/n(i+1));
        else
            x(i)=0;
        end
    end
else
    n=1:length(D);
    n2=2:length(D);
   	x=zeros(1,length(n2));

    for i=1:length(n2)
        if D(i+1)>eps && D(i)>eps
            x(i)=(log(D(i+1))-log(D(i)))/log((n(i))/n(i+1));
        else
            x(i)=0;
        end    
    end

end
end

