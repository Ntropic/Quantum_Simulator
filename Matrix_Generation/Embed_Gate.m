function [ A_new cluster_list ] = Embed_Gate( A ,indexes,N )
%EMBED_GATE embeds the transformation matrix A (with size(A)=2^len(indexes)*[1 1])
% of the qubits with indexes into an N qubit quantum computer


if isa(A,'sym')==0 %| arg_num==0
    
    if isa(A,'function_handle')==0 %is numeric - not a function handle
        
            sizle=size(A);
            leni=length(indexes);

            if sizle(1)~=sizle(2)
                error('A has to be a square matrix.')
            end
            if sizle(1)~=2^leni
                error('A has to be a matrix of size 2^(length(indexes)).')
            end
            if N<leni
                error('N has to be larger or equal to length(indexes).')
            end

        %Determine the datatype of A
        if isa(A,'sym')==1 
           arg_num=length(symvar(A));
        end
        if any(indexes>N)
            error('Index can not be higher than number of qbits.')
        end

        A_new=sparse(2^N,2^N);


        indexes=fliplr(indexes);
        other_indexes=setdiff(1:N,indexes);

        %Create blockstructure for the repetitions
        %1. the qubits of the indexes list shall be zero. the possibilities under
        %this boundary condition are the starting points for the blocks.
        %2. the possible numbers of the indexes lists qubits under
        %the boundary condition of all other qubits being zero are added to 
        %every starting point to create the block we then input the matrix A
        %into each block

        %1. Create block starting points by setting the qubits in indexes to zero
        perms=zeros(2^length(other_indexes),N);
        perms(:,other_indexes)=dec2bin(0:(2^length(other_indexes)-1),length(other_indexes))-'0';
        st_po=bi2de(perms)+1;

        %2. Find numbers, that should be added to starting points to get blocks
        perms=zeros(2^length(indexes),N);
        perms(:,indexes)=dec2bin(0:(2^length(indexes)-1),length(indexes))-'0';
        ad_po=bi2de(perms);

        cluster_list=repmat(st_po,1,length(ad_po))+repmat(ad_po',length(st_po),1);

        for i=1:length(st_po)
            A_new(cluster_list(i,:),cluster_list(i,:))=A;
        end
    else    %isa functionn_handle
        
        [ c_mat,c_vars ] = fun_mat2str( A );

            sizle=size(c_mat);
            leni=length(indexes);

            if sizle(1)~=sizle(2)
                error('A has to be a square matrix.')
            end
            if sizle(1)~=2^leni
                error('A has to be a matrix of size 2^(length(indexes)).')
            end
            if N<leni
                error('N has to be larger or equal to length(indexes).')
            end

        A_new=cell(2^N,2^N);
        A_new(:)={'0.0'}; 
        
        indexes=fliplr(indexes);
        other_indexes=setdiff(1:N,indexes);

        %Create blockstructure for the repetitions
        %1. the qubits of the indexes list shall be zero. the possibilities under
        %this boundary condition are the starting points for the blocks.
        %2. the possible numbers of the indexes lists qubits under
        %the boundary condition of all other qubits being zero are added to 
        %every starting point to create the block we then input the matrix A
        %into each block

        %1. Create block starting points by setting the qubits in indexes to zero
        perms=zeros(2^length(other_indexes),N);
        perms(:,other_indexes)=dec2bin(0:(2^length(other_indexes)-1),length(other_indexes))-'0';
        st_po=bi2de(perms)+1;

        %2. Find numbers, that should be added to starting points to get blocks
        perms=zeros(2^length(indexes),N);
        perms(:,indexes)=dec2bin(0:(2^length(indexes)-1),length(indexes))-'0';
        ad_po=bi2de(perms);

        cluster_list=repmat(st_po,1,length(ad_po))+repmat(ad_po',length(st_po),1);

        for i=1:length(st_po)
            A_new(cluster_list(i,:),cluster_list(i,:))=c_mat;
        end      
        [ A_new ] = str2fun_mat( A_new,c_vars );
        
    end
    
else        %symbolic
    
        sizle=size(A);
        leni=length(indexes);

        if sizle(1)~=sizle(2)
            error('A has to be a square matrix.')
        end
        if sizle(1)~=2^leni
            error('A has to be a matrix of size 2^(length(indexes)).')
        end
        if N<leni
            error('N has to be larger or equal to length(indexes).')
        end

        %Determine the datatype of A
        if isa(A,'sym')==1 
           arg_num=length(symvar(A));
        end
        if any(indexes>N)
            error('Index can not be higher than number of qbits.')
        end

    A_new=sym(zeros(2^N,2^N));
    indexes=fliplr(indexes);
    other_indexes=setdiff(1:N,indexes);

    %Create blockstructure for the repetitions
    %1. the qubits of the indexes list shall be zero. the possibilities under
    %this boundary condition are the starting points for the blocks.
    %2. the possible numbers of the indexes lists qubits under
    %the boundary condition of all other qubits being zero are added to 
    %every starting point to create the block we then input the matrix A
    %into each block

    %1. Create block starting points by setting the qubits in indexes to zero
    perms=zeros(2^length(other_indexes),N);
    perms(:,other_indexes)=dec2bin(0:(2^length(other_indexes)-1),length(other_indexes))-'0';
    st_po=bi2de(perms)+1;

    %2. Find numbers, that should be added to starting points to get blocks
    perms=zeros(2^length(indexes),N);
    perms(:,indexes)=dec2bin(0:(2^length(indexes)-1),length(indexes))-'0';
    ad_po=bi2de(perms);

    cluster_list=repmat(st_po,1,length(ad_po))+repmat(ad_po',length(st_po),1);

    for i=1:length(st_po)
        A_new(cluster_list(i,:),cluster_list(i,:))=A;
    end
    
end
end