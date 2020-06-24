function [ block ] = Block_Parallelizer( block )
%BLOCK_PARALLELIZER is used in creating circuit diagrams of quantum
%circuits. It is called by Gates2Tex.m. It receives a block and searches
%for gates that can be parallelized, parallelizes the block and returns a
%smaller block. 

B=1-strcmp(block,'\\qw');
C=zeros(size(B));

ind_list=zeros(2,size(B,2));

for i=1:size(B,2)
    f=find(B(:,i));
    ind_list(:,i)=[min(f);max(f)];
    C(min(f):max(f),i)=1;
end

s_or=size(block,2);
j=2;
while j<=size(block,2)
    if j<s_or %symmetric approach
        si=0;
        si1=0;
        if length(find(C(ind_list(1,j):ind_list(2,j),j-1)))==0 %From left to the right
            %How many can be found?
            how_many=1;
            truer=1;
            while truer && j-how_many-1>=1 && j+how_many<=size(C,2)
                truer1=isequal(C(:,j-how_many-1),C(:,j-1));
                truer2=isequal(C(:,j+how_many),C(:,j));
                truer=(truer1 && truer2);
                if truer
                    how_many=how_many+1;
                end
            end
            %Space for parallelization is there
            for k=1:how_many
                f=find(B(:,j+k-1));
                block(f,j-how_many+k-1)=block(f,j+k-1);
                C(f,j-how_many+k-1)=C(f,j+k-1);
                B(f,j-how_many+k-1)=B(f,j+k-1);
            end
            ind_list(:,j:(j+k-1))=[];
            block(:,j:(j+k-1))=[];
            C(:,j:(j+k-1))=[];
            B(:,j:(j+k-1))=[];
            for i=1:size(B,2)
                f=find(B(:,i));
                ind_list(:,i)=[min(f);max(f)];
                C(min(f):max(f),i)=1;
            end
        else
            si=1;
        end
        
        j2=size(B,2)-j+1;
        if j2>1
            if length(find(C(ind_list(1,j2):ind_list(2,j2),j2-1)))==0 %From left to the right
                %How many can be found?
                how_many=1;
                truer=1;
                while truer && j2-how_many-1>=1 && j2+how_many<=size(C,2)
                    truer1=isequal(C(:,j2-how_many-1),C(:,j2-1));
                    truer2=isequal(C(:,j2+how_many),C(:,j2));
                    truer=(truer1 && truer2);
                    if truer
                        how_many=how_many+1;
                    end
                end
                %Space for parallelization is there
                for k=1:how_many
                    f=find(B(:,j2+k-1));
                    block(f,j2-how_many+k-1)=block(f,j2+k-1);
                    C(f,j2-how_many+k-1)=C(f,j2+k-1);
                    B(f,j2-how_many+k-1)=B(f,j2+k-1);
                end
                ind_list(:,j2:(j2+k-1))=[];
                block(:,j2:(j2+k-1))=[];
                C(:,j2:(j2+k-1))=[];
                B(:,j2:(j2+k-1))=[];
                for i=1:size(B,2)
                    f=find(B(:,i));
                    ind_list(:,i)=[min(f);max(f)];
                    C(min(f):max(f),i)=1;
                end
            else
                si1=1;
            end
        else
            si1=1;
        end

        if si==1 && si1==1
            j=j+1;
        end
    else
%         if length(find(C(ind_list(1,j):ind_list(2,j),j-1)))==0
%             %Space for parallelization is there
%             f=find(B(:,j));
%             block(f,j-1)=block(f,j);
%             fprintf('c')
%             C(f,j-1)=C(f,j)
%             B(f,j-1)=B(f,j)
%             ind_list(:,j)=[];
%             block(:,j)=[];
%             C(:,j)=[];
%             B(:,j)=[];
%         else
%             j=j+1;
%         end
        j=j+1;
    end
end
end

