%draw_fractal_bintree_bin_code_list.m
clc;
clear all;
close all;
fclose('all');

p=pwd;
elem=strsplit(p,'\');
shortened=strjoin(elem((1:length(elem)-3)),'\');
addpath(genpath(shortened));

%Create binary tree
depth_max=7;
n=1;

f=figure();
for depth=6%2:depth_max
    tree=bintree(depth,n);

    perm=all_num_algorithm(tree);
    perm=perm-1;
    bin_perm=(dec2bin(perm)-'0')';
    bin_perm=bin_perm(end:-1:1,:);
    bin_perm=[bin_perm;zeros(1,size(bin_perm,2))];
    bin_perm=[bin_perm zeros(size(bin_perm,1),1)];
    
    %Plot
    ax1=subplot(2,1,1);
    x0=repmat((1:size(perm,2))-1,2,1);
    x0(1,:)=x0(1,:)-0.5;
    x0(2,:)=x0(2,:)+0.5;
    x=x0(:)';
    p=repmat(perm,2,1);
    plot(x,p(:)','k');
    axis tight;
    
    ax2=subplot(2,1,2);
    pcolor(1-bin_perm);
    g=gray(5);
    colormap(g([2,5],:));
    ax2.YTick=(1:depth)+0.5;
    a=ax2.YTickLabel;
    for i=1:length(a)
        a{i}=a{i}(1);
    end
    ax2.YTickLabel=a;
    s=ceil(size(perm,2)/15);
    ax2.XTick=(1:s:size(perm,2))+0.5;
    a={};
    m=1;
    for i=1:s:size(perm,2)
        a{m}=num2str(i-1);
        m=m+1;
    end
    ax2.XTickLabel=a;
    ax1.XTick=0:s:size(perm,2)-1;
    ax1.XTickLabel=ax2.XTickLabel;
    pause(1)
end

matlab2tikz('perm_bin_6.tex')