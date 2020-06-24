function [ block ] = Block_Corrector( block,indexes,blocker2,indexer2,circuit_subs )
%BLOCK_CORRECTOR Corrects block elements 
    % Check for multigates and deal with them accordingly
    if nargin>2
        if length(blocker2)>0
            block2=blocker2;
            indexes2=indexer2;
        end
    end
    stringer='multigate';
    for i=1:size(block,2) 
        pos=0;
        for j=1:size(block,1)
            %Find the multigate marker
            p=strfind(block{j,i},stringer);
            if length(p)>0
                pos=p;
                posj=j;
            end
        end
        if pos>0    %Found multigate
            b=block{posj,i};
            stringer2=b(1:end-1);
            stringer2=strsplit(stringer2,'}{');
            stringer2=stringer2(2:end);
            stringer2=strjoin(stringer2,'}{');
            
            %Find other gates of this type
            ind_list=[];
            for j=1:size(block,1)
                p=strfind(block{j,i},stringer2);
                if length(p)>0
                   ind_list=[ind_list j];
                end
            end
            %Check indexes (are they together?)
            in_cmp=min(ind_list):max(ind_list);
            if length(ind_list)==length(in_cmp)  %multigate is together
                %Reorder multigate elements (to get a compiling code)
                [block{in_cmp(1:end-1),i}]=deal(['\\ghost{' stringer2 '}']);
                block{in_cmp(end),i}=['\\multigate{' num2str(length(in_cmp)-1) '}{' stringer2 '}'];
            else % multigate is split -> take alternative representation
                if exist('block2')
                    %Exchange elements with
                    b=block2{posj,i};
                    if strfind(b,'multigate');
                        stringer2=b(1:end-1);
                        stringer2=strsplit(stringer2,'}{');
                        stringer2=stringer2(2:end);
                        stringer2=strjoin(stringer2,'}{');
                    else
                        stringer2=b(1:end-1);
                        stringer2=strsplit(stringer2,'{');
                        stringer2=stringer2(2:end);
                        stringer2=strjoin(stringer2,'{');
                    end
                    in_cmp=sort(ind_list,2,'descend');
                    diff_in_cmp=-diff(in_cmp);
%% Old approach using sgates
%                     for j=length(in_cmp):-1:2
%                          block{in_cmp(j),i}=['\\sgate{' stringer2 '}{' num2str(in_cmp(j)-in_cmp(j-1)) '}'];
%                     end
%                     block{in_cmp(1),i}=['\\gate{' stringer2 '}'];

%% New approach (using multigates)
                    if diff_in_cmp(1)>1
                        block{in_cmp(1),i}=['\\gate{' stringer2 '}\\qwx[' num2str(diff_in_cmp(1)) ']'];
                    else
                        %Count the number of repeating ones
                        o=1;
                        while diff_in_cmp(o+1)==1
                            o=o+1;
                        end
                        block{in_cmp(1),i}=['\\multigate{' num2str(o) '}{' stringer2 '}'];
                    end
                    for j=2:length(in_cmp)
                        %Check if belongs to a multigate
                        if diff_in_cmp(j-1)==1 %part of a multigate
                            if j<length(in_cmp)
                                if diff_in_cmp(j)==1 %Multigate continues
                                    block{in_cmp(j),i}=['\\ghost{' stringer2 '}'];
                                else %Multigate ends
                                    block{in_cmp(j),i}=['\\ghost{' stringer2 '}\\qwx[' num2str(diff_in_cmp(j)) ']'];
                                end
                            else
                                block{in_cmp(j),i}=['\\ghost{' stringer2 '}'];
                            end
                        else %Start new block
                            if j<length(in_cmp)
                                if diff_in_cmp(j)==1 %Multigate
                                    %Count the number of repeating ones
                                    o=1;
                                    if j+o<length(diff_in_cmp)
                                        while diff_in_cmp(j+o+1)==1 && j+o<length(diff_in_cmp)
                                            o=o+1;
                                        end
                                    end
                                    block{in_cmp(j),i}=['\\multigate{' num2str(o) '}{' stringer2 '}'];
                                else %Not a Multigate
                                    block{in_cmp(j),i}=['\\gate{' stringer2 '}\\qwx[' num2str(diff_in_cmp(j)) ']'];
                                end
                            else
                                block{in_cmp(j),i}=['\\gate{' stringer2 '}'];
                            end
                        end
                    end
                    
                    
                    
                    
                else
                    %in_cmp=sort(ind_list)
                    in_cmp=sort(ind_list,2,'descend');
                    diff_in_cmp=-diff(in_cmp);
%% Old approach using sgates
%                     for j=length(in_cmp):-1:2
%                          block{in_cmp(j),i}=['\\sgate{' stringer2 '}{' num2str(in_cmp(j)-in_cmp(j-1)) '}'];
%                     end
%                     block{in_cmp(1),i}=['\\gate{' stringer2 '}'];

%% New approach (using multigates)
                    if diff_in_cmp(1)>1
                        block{in_cmp(1),i}=['\\gate{' stringer2 '}\\qwx[' num2str(diff_in_cmp(1)) ']'];
                    else
                        %Count the number of repeating ones
                        o=1;
                        while diff_in_cmp(o+1)==1
                            o=o+1;
                        end
                        block{in_cmp(1),i}=['\\multigate{' num2str(o) '}{' stringer2 '}'];
                    end
                    for j=2:length(in_cmp)
                        %Check if belongs to a multigate
                        if diff_in_cmp(j-1)==1 %part of a multigate
                            if j<length(in_cmp)
                                if diff_in_cmp(j)==1 %Multigate continues
                                    block{in_cmp(j),i}=['\\ghost{' stringer2 '}'];
                                else %Multigate ends
                                    block{in_cmp(j),i}=['\\ghost{' stringer2 '}\\qwx[' num2str(diff_in_cmp(j)) ']'];
                                end
                            else
                                block{in_cmp(j),i}=['\\ghost{' stringer2 '}'];
                            end
                        else %Start new block
                            if j<length(in_cmp)
                                if diff_in_cmp(j)==1 %Multigate
                                    %Count the number of repeating ones
                                    o=1;
                                    if j+o<length(diff_in_cmp)
                                        while diff_in_cmp(j+o+1)==1 && j+o<length(diff_in_cmp)
                                            o=o+1;
                                        end
                                    end
                                    block{in_cmp(j),i}=['\\multigate{' num2str(o) '}{' stringer2 '}'];
                                else %Not a Multigate
                                    block{in_cmp(j),i}=['\\gate{' stringer2 '}\\qwx[' num2str(diff_in_cmp(j)) ']'];
                                end
                            else
                                block{in_cmp(j),i}=['\\gate{' stringer2 '}'];
                            end
                        end
                    end    
                    
                end
            end
        end
    end
    
    
    %Exchange placeholder elements by appropriate elements
    for i=1:size(indexes,2)
        %Search for _i_ and replace it with the difference j-i
        stringer=['_' num2str(i) '_'];
        ind_i=indexes(1,i);
        for j=ind_i+1:size(block,1)
            %Does index belong to a multiblock? which element is the closest one?
            if exist('ind_list')
                %Connect to closest element in multiblock
                an=any(ind_list==ind_i);
                if an==1
                    [x y]=min(abs(j-ind_list));
                    ind_i=ind_list(y);
                end
            end
            replacer=num2str(j-ind_i);

            for k=1:size(block,2)
                if length(strfind(block{j,k},stringer))>0
                    ind_i=j;
                end
            end
            block(j,:)=cellfun(@(x) strrep(x,stringer,replacer),block(j,:),'UniformOutput',false);
        end
        stringer=['_' num2str(i) '_'];
        ind_i=indexes(1,i);
        for j=ind_i:-1:1
            %Does index belong to a multiblock? which element is the closest one?
            if exist('ind_list')
                %Connect to closest element in multiblock
                an=any(ind_list==ind_i);
                if an==1
                    [x y]=min(abs(j-ind_list));
                    ind_i=ind_list(y);
                end
            end
            replacer=num2str(j-ind_i);

            for k=1:size(block,2)
                if length(strfind(block{j,k},stringer))>0
                    ind_i=j;
                end
            end
            block(j,:)=cellfun(@(x) strrep(x,stringer,replacer),block(j,:),'UniformOutput',false);
        end
    end
    
    %Circuit_subs
    if length(circuit_subs)>0
        circuit_subs=strrep(circuit_subs,'\','\\');
        for i=1:length(circuit_subs)/2
            for j=1:size(block,1)
                block(j,:)=cellfun(@(x) strrep(x,circuit_subs{2*i-1},circuit_subs{2*i}),block(j,:),'UniformOutput',false);
            end
        end
    end
end

