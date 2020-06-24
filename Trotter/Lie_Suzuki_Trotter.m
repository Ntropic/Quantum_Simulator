function [ T,j_list ] = Lie_Suzuki_Trotter( Uij,l,print,method,permutations )
%LIE_SUZUKI_TROTTER does the Lie-Suzuki-Trotter decomposition 
vec_size=size(Uij{1},1);
if nargin==3
    method='ordered';
elseif nargin==2
    print='y'
    method='ordered';
end



if print=='y'
    fprintf('->Do Lie-Suzuki-Trotter Steps\n');
end
if length(strfind(method,'ordered'))>0
    j_list=1:length(Uij);
    T=Uij{1};
    for j=2:length(Uij)
        T=T*Uij{j};
    end
    T=T^l;
    j_list=repmat(j_list,1,l);
elseif length(strfind(method,'random'))>0
    if length(strfind(method,'every'))==0
        j_list=randperm(length(Uij));
        T=Uij{j_list(1)};

        for j=2:length(Uij)
            T=T*Uij{j_list(j)};
        end
        T=T^l;
        j_list=repmat(j_list,1,l);
    else
        j_list=[];
        for i=1:l
            j_list_new=randperm(length(Uij));
            if i==1
                T=Uij{j_list_new(1)};
            end
            
            for j=2:length(Uij)
                T=T*Uij{j_list_new(j)};
            end
            j_list=[j_list j_list_new];
        end
    end
elseif length(strfind(method,'permutations'))>0
    if size(permutations,1)==1
        if length(permutations)==length(Uij)
            T=Uij{permutations(1)};

            for j=2:length(Uij)
                T=T*Uij{permutations(j)};
            end
            T=T^l;
            j_list=repmat(permutations,1,l);
        elseif length(permutations)>length(Uij)
            T=Uij{permutations(1)};
            j_list=permutations;
            j_list=repmat(permutations,1,ceil(l*length(Uij)/length(permutations)));
            j_list=j_list(1:length(Uij)*l);

            for j=2:length(Uij)*l
                T=T*Uij{permutations(j)};
            end
        end
    elseif size(permutations,1)>1
        if length(permutations)==length(Uij)
            T=Uij{permutations(1)};

            for j=2:length(Uij)
                T=T*Uij{permutations(j)};
            end
            T=T^l;
            j_list=repmat(permutations,1,l);
        elseif length(permutations)>length(Uij)
            T=Uij{permutations(1)};
            j_list=permutations;
            j_list=repmat(permutations,1,ceil(l*length(Uij)/length(permutations)));
            j_list=j_list(1:length(Uij)*l);

            for j=2:length(Hij)*l
                T=T*Uij{permutations(j)};
            end
        end
    end
end
end

