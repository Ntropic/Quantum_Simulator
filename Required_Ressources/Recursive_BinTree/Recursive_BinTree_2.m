function [ perm ] = Recursive_BinTree_2( N,m,i,perm )
%RECURSIVE_BINTREE creates a permutation of the natural numbers from m to 
%2^N-1, where subsequent numbers have a Hamming distance of 1
% N     -> defines 2^N-1 (upper limit)
% m     -> defines lower limit
% i     -> current level (bit on which to act from right to left)
% perm  -> sequence (up to this point)


    if i>1  %Not on the lowest level
        perm=Recursive_BinTree_2(N,m,i-1,perm);     %Call function on lower level
    end
    p_end=perm(end);                        %Last element of sequence
    bin_i=mod(floor(p_end/(2^(i-1))),2);      %Intersection bit of last 
                                              %element in sequence
    p_end=p_end+2^(i-1)*(1-2*bin_i);          %Flip i'th bit
    if p_end>m                                %Is new element is in range?
        perm=[perm p_end];                    %Add new (flipped) element
        if i>1  %Not on lowest level
            perm=Recursive_BinTree_2(N,m,i-1,perm); %Call function on lower 
                                                  %level
        end
    end
end

