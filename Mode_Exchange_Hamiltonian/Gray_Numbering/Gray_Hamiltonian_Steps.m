function [ Hij ] = Gray_Hamiltonian_Steps( n,mode )
%Creates a list of Hamiltonians H for the xx+yy interactions of a Gray coded 
%Boson Mode Interaction
if nargin==1
    mode=0;
elseif nargin==2
%     if length(strfind(mode,'ladder'))>0 %All expansions are now using
%     even uneven ordering
%         mode=1;
%     else
%         mode=0;
%     end
end

% if mode==0 %Typical Gray ordering
    %% Prepare elements
    s=ceil(log2(n+1));

    tree=bintree(s,1);
    gray=all_num_algorithm(tree);

    gray_ind_list=[];
    prefactor=[];

    for i=[1:2:n 2:2:n]
        for l=1:n
            j=l-i+1;
            if j>=1
                %prefactor=[prefactor;1+0*(i*j)];  %Strength of interaction 
                prefactor=[prefactor;sqrt(i*j)];

                gray_num=gray([i+1,j]);
                graybin=dec2bin(gray_num-1,s)';
                graybin=graybin(:)';
                graybin=graybin-'0';

                gray_num2=gray([i,j+1]);
                graybin2=dec2bin(gray_num2-1,s)';
                graybin2=graybin2(:)';
                graybin2=graybin2-'0';
                gray_ind_list=[gray_ind_list; bin2dec(num2str(graybin))+1,bin2dec(num2str(graybin2))+1];
            end
        end
     end   

    %% Create the Hamiltonians
    H=sparse(2^(2*s),2^(2*s));
    for o=1:size(gray_ind_list,1)
        H_step=H;
        H_step(gray_ind_list(o,1),gray_ind_list(o,2))=prefactor(o);
        H_step(gray_ind_list(o,2),gray_ind_list(o,1))=prefactor(o);
        Hij{o}=H_step;
    end
    
    %Order even and uneven terms
    Hij={Hij{1:2:length(Hij)},Hij{2:2:length(Hij)}};
% 
% elseif mode==1 %% Order operations to Ladder ordering, but keeping the gray coding to save qubits
%     %% Prepare interactions
%     s=ceil(log2(n+1));
% 
%     tree=bintree(s,1);
%     gray=all_num_algorithm(tree);
% 
%     block_nums=[];
%     prefactor=[];
%     gray_ind_list=[];
%     for d=0:n
%         i_max=floor((n+d+1)/2);
%         for i=d+1:i_max
%             j=i-d;
%             block_nums=[block_nums;i j];
%             prefactor=[prefactor;sqrt(i*j)];  %Strength of interaction 
%         end
%         if d~=0
%             for i=d+1:i_max
%                 j=i-d;
%                 block_nums=[block_nums;j i];
%                 prefactor=[prefactor;sqrt(i*j)];  %Strength of interaction 
%             end
%         end
%     end
%     for n_now=1:n
%         for i=1:2:n_now
%             j=n_now+1-i;
%             block_nums=[block_nums;i j];
%             prefactor=[prefactor;sqrt(i*j)];  %Strength of interaction 
%         end
%         for i=2:2:n_now
%             j=n_now+1-i;
%             block_nums=[block_nums;i j];
%             prefactor=[prefactor;sqrt(i*j)];  %Strength of interaction 
%         end
%     end
% 
%     for i=1:size(block_nums,1)
%         gi=2^s*(gray(block_nums(i,1))-1)+gray(block_nums(i,2)+1);
%         gf=2^s*(gray(block_nums(i,1)+1)-1)+gray(block_nums(i,2));
%         
%         gray_ind_list=[gray_ind_list; gi,gf];
%     end
% 
%     %% Create the Hamiltonians
%     H=sparse(2^(2*s),2^(2*s));
%     for o=1:size(block_nums,1)
%         H_step=H;
%         H_step(gray_ind_list(o,1),gray_ind_list(o,2))=prefactor(o);
%         H_step(gray_ind_list(o,2),gray_ind_list(o,1))=prefactor(o);
%         Hij{o}=H_step;
%     end
end

