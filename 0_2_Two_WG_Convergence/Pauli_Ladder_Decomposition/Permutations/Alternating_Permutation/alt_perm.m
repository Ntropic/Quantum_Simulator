function [ Am,pm ] = alt_perm( m_max )
%Calculates the number of permutations of length m accoding to the rules of
%not changing nearest neighbor commutation

Am=1; %Initialize
pm=1; %Initialize
for m=2:m_max
    pm_1=pm;
    pm=zeros(1,m); 
    for j=1:m-1
        pm(j)=sum(pm_1(1:m-j)); %Determine entringer numbers
    end
    Am(m)=sum(pm); %Determine Euler Up-Down numbers
end

end

