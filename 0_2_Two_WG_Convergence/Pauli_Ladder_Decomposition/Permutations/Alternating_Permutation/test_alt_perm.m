%test_alt_perm.m
clc;
clear all;
close all;

m_max=15;

mode=3;

if mode==1
    %for m=1:m_max
        [Am,pm]=alt_perm(m_max)
    %    Am
    %    pm
    %    fprintf('--------------\n')
    %end
elseif mode==2
    Am=alt_perm(m_max);  
    plot(1:m_max,Am./factorial(1:m_max),'r');
    %hold on;
    %plot(1:m_max,factorial(1:m_max),'b');
else
    [Am,pm]=alt_perm(m_max);
    n=2:m_max;
    2*n.*Am(n-1)./Am(n)

end
