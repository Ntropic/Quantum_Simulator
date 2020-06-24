function [ indexes ] = Rec_Multi_Gray_Indexes( m,n,s )
%REC_MULTI_LADDER_INDEXES creates a list of indexes in the subspace of the 
%ladder coding for m modes it is used within Multi_Ladder_Indexes

ladder_red=2.^(0:n);

s_m=2^(s*(m-1));
indexes=[];
if m>2
    for i=1:n+1
        indexes=[indexes,s_m*ladder_red(i)+Rec_Multi_Ladder_Indexes( m-1,n-i+1,s)];
    end
else
    for i=1:n+1
        for j=1:n+2-i
            indexes=[indexes,s_m*ladder_red(i)+ladder_red(j)];
        end
    end
end
end

