function [ A,B ] = testrandomsymgenerator( symlist,dim1,dim2,prob )
x=symlist(1);
y=symlist(2);
c=symlist(3);
d=symlist(4);

elements={x,x^2,cos(x),sin(x),sin(y),cos(y),y,y^2,c,d};

A=sym(zeros(dim1));
B=sym(zeros(dim2));

for i=1:dim1(1)
    for j=1:dim1(2)
        if rand(1)>1-prob
            A(i,j)=elements{ceil(length(elements)*rand(1))};
        end
    end
end
for i=1:dim1(1)
    for j=1:dim1(2)
        if rand(1)>1-prob
            B(i,j)=elements{ceil(length(elements)*rand(1))};
        end
    end
end

end

