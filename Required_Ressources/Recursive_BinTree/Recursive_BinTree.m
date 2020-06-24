function [ perm ] = Recursive_BinTree( N0,N,m,varargin )
%RECURSIVE_BINTREE creates a permutation of the natural numbers m to 2^M
%that has a Hamming distance of 1 between adjacent elements within the
%permutation
if nargin==4
    m_bin=varargin{1};
    if N==1
        b=m_bin(length(m_bin)-N+1);
        perm=[m_bin;m_bin];
        perm(2,end)=mod(b+1,2);
    else
        perm=Recursive_BinTree(N0,N-1,m,m_bin);
        p_end=perm(end,:);
        p_end(N0-N+1)=mod(p_end(N0-N+1)+1,2);
        perm2=Recursive_BinTree(N0,N-1,m,p_end);
        perm=[perm;perm2];
    end
elseif nargin==3
    m_bin=dec2bin(m-1,N0)-'0';
    if N>1
        m_bin=m_bin(N0-N+1);
        perm=Recursive_BinTree(N0,N-1,m);
        if m_bin==0
            p_end=perm(end,:);
            p_end(N0-N+1)=mod(p_end(N0-N+1)+1,2);
            perm2=Recursive_BinTree(N0,N-1,m,p_end);
            perm=[perm; perm2];
        end
    else
        perm=m_bin;
        if m_bin(end)==0
            perm2=perm;
            perm2(end)=1;
            perm=[perm; perm2];
        end
    end
end
end

