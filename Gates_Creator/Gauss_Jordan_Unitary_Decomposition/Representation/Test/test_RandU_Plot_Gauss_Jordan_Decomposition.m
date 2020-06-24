%test_Plot_Gauss_Jordan_Decomposition.m
clear all;
close all;
clc;

make_mov=1;

p=pwd;
if any(strfind(p,'\'));
    elem=strsplit(p,'\');
else
    elem=strsplit(p,'/');
end
shortened=fullfile(elem{1:end-4});
addpath(genpath(shortened));

%Number of qubits per mode
s=3;
file_name=['randu_gj_' num2str(s) ];
matrix=RandU(2^s);


%% Default
mat=Plot_Gauss_Jordan_Decomposition(matrix,1); 

h=figure('Position',[100 100 600 600]);
ax1=subplot(1,1,1);

for i=1:length(mat)
    name_now=['r_gray_none_tex_' num2str(i)];
    [bra_str ket_str indexes]=Add_PColorMat(ax1,full(mat{i}) ,{name_now,s,2});
    
    drawnow();
    B=PColor2Tikz2Pdf2Png(bra_str,ket_str,indexes,'PColor_Tikz',name_now,600);
    
end