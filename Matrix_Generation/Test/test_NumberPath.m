%test_NumberPath.m
clc;
clear all;
close all;

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-2)),'\');
addpath(genpath(shortened))

m=8;
%Initialize
index_list=NumberPaths(m);
surf(index_list,'EdgeColor','none')
view(0,90)
%Calculate Path

%Plot Results
%phi=linspace(0,pi,50);
%max_r=2^3;
%cmap=jet(max_r*2);
%for i=1:length(index_list)
%    b=index_list(i,1);
%    m=mean(b);
%    r=abs((b(2)-(1))/2);
%    x=m+r*cos(phi);
%    y=abs(r*sin(phi));
%    plot(x(1:end-1),y(1:end-1),'Color',cmap(r*2,:));
%    quiver(x(end-1),y(end-1),x(end)-x(end-1),y(end)-y(end-1),0, 'MaxHeadSize',0.5,'Color',cmap(r*2,:));
%    hold on;
%end
%axis equal;
%axis([0 max(index_list(:))+1 0 (max(index_list(:)))/4+0.5]);