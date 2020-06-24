function [ h ] = permutation_draw( perm,mode,h )
%PERMUTATION_DRAW draws the permutation list (perm) according to mode

if nargin==1
    mode='circ';
    h=figure();
    a=axes(h);  
    %set(a,'Visible','off');
elseif nargin==2
    h=figure();
    a=axes(h);
    %set(a,'Visible','off');
elseif nargin==3
    a=get(h,'children');
    a=a(1);
    %set(a,'Visible','off');
end

perm=perm(:);
if length(strfind(mode,'circ'))>0
    phi=linspace(0,pi,50);
    cmap=jet(length(perm)-1);
    
    perm_list=[perm(1:end-1),perm(2:end)];
    for i=1:length(perm)-1
        b=perm_list(i,:);
        m=mean(b);
        r=(max(b)-min(b))/2;
        x=m+r*cos(phi);
        y=r*sin(phi);
        plot(a,x,y,'Color',cmap(i,:));
        hold on;
    end
    axis equal;
    axis([min(perm) max(perm) 0 max(abs(diff(perm)))/2*1.1]);
elseif length(strfind(mode,'lines1'))>0
    y1=-(1:size(perm,1));
    x1=perm;
    for i=1:length(x1)-1
        x=x1(i:i+1);
        y=y1(i)*[1 1];
        plot(a,x,y,'k');
        hold on;
    end
    axis equal;
    axis([1 max(perm) -length(perm)+1 1]);
elseif length(strfind(mode,'lines'))>0
    y=-(0:size(perm,1)-1);
    x=perm;
    plot(a,x,y,'k');
    axis equal;
    axis([1 max(perm) -length(perm)+1 0]);
end
end