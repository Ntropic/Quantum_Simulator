function [ F,bra ] = DensityBar3( rho,varargin )
%DENSITYBAR3 Creates a bar3 plot of a density matrix rho
%Optional:
%   Coloring    [bool] true,false
%   Fig         {f,ax} 
%   Bra         {'bra1','bra2'}
%   Ket         {'ket1','ket2'}
%   BraKetMode  {'fock_bin',photons,modes}
%   AbsValue    [bool] true,false
%   ZLim        [min,max]
%   ZLimTight   [bool] true,false

max_ticks=8;
tick_num=min([length(rho),max_ticks]);
div=length(rho)/tick_num;  
indexes=1:div:length(rho);

if nargin>1
    var_len=length(varargin);
    if mod(var_len,2)~=0
        error('Variable length input must contain variables for every argument!');
    end
    var={};
    arg={};
    for i=1:var_len
        if mod(i,2)==0
            var={var{:},varargin{i}};
        else
            arg={arg{:},varargin{i}};
        end
    end
else 
    arg={};
end

%% Determine auxiliary variables
Coloring=true;
Fig=false;
bra_given=false;
ket_given=false;
BraKetMode=false;
AbsValue=false;
z_lim_given=false;
z_lim_tight=false;

s=log2(size(rho,1));
mode={'fock_bin',1,s};

%% Check if these options are given within varargin
if length(find(strcmp(arg,'Coloring')))>0
    i=find(strcmp(arg,'Coloring'));
    Coloring=var{i};
end
if length(find(strcmp(arg,'Fig')))>0
    i=find(strcmp(arg,'Fig'));
    Fig=true;
    f=var{i}{1};
    ax=var{i}{2};
end
if length(find(strcmp(arg,'Bra')))>0
    i=find(strcmp(arg,'Bra'));
    bra_given=true;
    bra=var{i};
end
if length(find(strcmp(arg,'Ket')))>0
    i=find(strcmp(arg,'Ket'));
    ket_given=true;
    ket=var{i};
end
if length(find(strcmp(arg,'BraKetMode')))>0
    i=find(strcmp(arg,'BraKetMode'));
    BraKetMode=true;
    mode=var{i};
end
if BraKetMode==true | bra_given==false | ket_given==false
    bra=bra_fock_tex(1:length(rho),mode);
    ket=ket_fock_tex(1:length(rho),mode);
    for i=1:length(bra)
        bra{i}=['$',bra{i}(2:end),'$'];
        ket{i}=['$',ket{i}(2:end),'$'];
    end
end    
if length(find(strcmp(arg,'AbsValue')))>0
    i=find(strcmp(arg,'AbsValue'));
    AbsValue=var{i};
end   
if length(find(strcmp(arg,'ZLim')))>0
    i=find(strcmp(arg,'ZLim'));
    z_lim_given=true;
    z_lim=var{i};
end   
if length(find(strcmp(arg,'ZLimTight')))>0
    i=find(strcmp(arg,'ZLimTight'));
    z_lim_tightval{i}=var{i};
end   

%% Plot
if Fig==false
    f=figure();
    ax=axes(f);
end
if max(abs(imag(rho(:))))==0
    if AbsValue==true;
        b=bar3(abs(rho(end:-1:1,:)),1);
    else
        b=bar3(rho(end:-1:1,:),1);
    end
else
    if AbsValue==true;
        b=bar3(abs(rho(end:-1:1,:)),1);
    else
        b=bar3(real(rho(end:-1:1,:)),1);
    end
end

if Coloring==true
    for k = 1:length(b)
        zdata = b(k).ZData;
        b(k).CData = zdata;
        b(k).FaceColor = 'interp';
    end
end



zlabel(ax,'p');
ax.XTick=indexes;
ax.XTickLabel=ket(indexes);

ax.YTick=indexes;
ax.YTickLabel=bra(indexes(end:-1:1));


if z_lim_given==false
    if AbsValue==true & z_lim_tight==false
        z=[0 1];
    elseif AbsValue==true & z_lim_tight==true
        z=[0 max(abs(rho(:)))];
    elseif AbsValue==false & z_lim_tight==true
        z=[min(abs(rho(:))*sign(rho(:))) max(abs(rho(:))*sign(rho(:)))];
    elseif AbsValue==false & z_lim_tight==false
        z=zlim;
    end
else
    z=z_lim;
end

caxis(z)

ax.XLim=[0.5 0.5+length(rho)];
ax.YLim=[0.5 0.5+length(rho)];
ax.ZLim=z;
asp=daspect;
asp(2)=asp(1);
daspect(asp);
view(45,45)

drawnow()
F={f,ax};
end

