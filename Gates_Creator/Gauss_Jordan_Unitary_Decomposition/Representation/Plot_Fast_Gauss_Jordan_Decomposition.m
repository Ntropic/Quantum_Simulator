function [ mat ] = Plot_Fast_Gauss_Jordan_Decomposition( matrix,plotting )
%GAUSS_JORDAN_DECOMPOSITION decomposes a a unitary matrix into a set of
%controlled U2 rotations -> in order to implement on a quantum computer
mini=10^-12;
n=log2(size(matrix,1));
    
if n==1         
    mat{1}=matrix;
    mat{2}=diag(ones(1,2^n));
else
    par1=sym('alpha');
    par2=sym('beta');
    par3=sym('delta');
    par4=sym('theta');
    par5=sym('phi');
    u_fun=@(alpha,beta,delta,theta)reshape([cos(theta/2)*exp(-(alpha*1i)/2-(beta*1i)/2+delta*1i) , sin(theta/2)*exp((alpha*1i)/2-(beta*1i)/2+delta*1i) , -sin(theta/2)*exp(-(alpha*1i)/2+(beta*1i)/2+delta*1i) , cos(theta/2)*exp((alpha*1i)/2 + (beta*1i)/2 + delta*1i) ],[2,2]);
    one=sparse(diag(ones(1,2^n)));
    c_bits=[];
    %Determine angles first, then reverse order and create gates
    angles_list=[];
    
    indexes2=[];
    bi=[];
    counter=1;
    for i=1:2^n           %Columns    
        %Find non zero elements in row i
        mat_i=matrix(:,i);
        ind=find(abs(mat_i)>mini);
        to_do=unique([i;ind(ind>=i)]);
        bits=dec2bin(to_do-1,n)-'0';
        intermed=bits(1,:);
        last=bits(1,:);
        for k=2:size(bits,1)
            differ=bits(k,:)-bits(k-1,:);
            for l=1:length(differ)
                if differ(l)~=0                
                    last(l)=last(l)+differ(l);
                    intermed=[intermed;last];
                end
            end
        end
        indexes=bin2dec(num2str(intermed,'%i'))+1;
                    
        how_many=0;
        if length(indexes)>1 && i<=2^n-1
            indexes2=[indexes2;indexes];
            for j=length(indexes):-1:2   %Rows
                index_now=[indexes(j),indexes(j-1)];
                %find pu index
                e=find(abs(intermed(j,:)-intermed(j-1,:)));
                indi=intermed(j,:);
                indi(e)=[];
                b=bin2dec(num2str(indi,'%i'))+1;
                pu_index_list=[b,e];
                %Create a zero on the element
                ind2zero=index_now(1);
                a=matrix(sort(index_now),i);

                if max(abs(a))>mini
                    u=sqrt(abs(a(1))^2+abs(a(2))^2);
                    U1=[a(2)/u,-a(1)/u;a(1)'/u,a(2)'/u];

                    if ind2zero~=min(index_now)
                        U1(1:2,:)=U1([2,1],:);
                        U1(2,:)=-U1(2,:);
                    end
                    
                    [alpha,beta,delta,theta]=Unitary2Angles(U1);
                    [alpha2,beta2,delta2,theta2]=Unitary2Angles(U1');
                    
                    %Check if angles are trivial
                    if max([alpha,beta,delta,theta])>mini
                        bi=[bi;ind2zero];
                        how_many=how_many+1;
                        d=dec2bin(indexes(j)-1,n)-'0';
                        c=d;
                        c(e)=NaN(1);   
                        c_bits=[c_bits;c];
                        
                        u=[alpha,beta,delta,theta];
                        u2=[alpha2,beta2,delta2,theta2];
                        pu_matrix=one;
                        pu_matrix(sort(index_now),sort(index_now))=u_fun(u(1),u(2),u(3),u(4));
                        matrix=pu_matrix*matrix;
                        %FockPrint(matrix)
                        mat{counter}=pu_matrix';
                        counter=counter+1;
                    end
                end
            end
        end
        if how_many==0
            s=angle(matrix(i,i));
            if abs(s)>mini
                %Rotate via controlled rotation
                pu_matrix=one;
                pu_matrix(i,i)=exp(-1i*s);
                matrix=pu_matrix*matrix;
                %FockPrint(matrix)   
                mat{counter}=pu_matrix';
                counter=counter+1;
            end
        end
    end
    mat{counter}=one;
end
%FockPrint(matrix)
%% Output 
mat=mat(end:-1:1);
if plotting==1    %Return cummulative matrix steps
    for i=2:length(mat)
        mat{i}=mat{i}*mat{i-1};
    end
end
end

