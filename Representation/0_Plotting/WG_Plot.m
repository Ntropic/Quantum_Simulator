function [ surf_mat x ti dx surf_mat_real surf_mat_imag] = WG_Plot( phi,t,inter,t_res,type )
%WG_Plot creates a surface plot of a WG system drawing the time evolution
%of the Waveguides
%   names - names of the basis states (remove unwanted components in
%           advance)
%   phi   - matrix [phi0,phi1,phi2,...,phin] of the wavefunctions at times 
%   t     - vector of times [t0,t1,...,tn]
%   inter - 0-non interpolated, 1-interpolated
%   t_res - resolution of the smallest intervall
%   type  - 0 -> sin , 1 -> exp & sin

if nargin==2
    inter=0;
    t_res=1;
    type=1;
end
if nargin==3
    t_res=1;
    type=1;
end

if size(phi,2)~=length(t)
    error('Every wavefunction needs a corresponding time')
end

%% Create vectors of the spatial distribution of a photon in a mode
len=size(phi,1);
res=50; %resolution per waveguide
dx=0.3; %width of waveguides <0.5
k=6;    %Decay factor

x=0.5:1/res:len+0.5; %position vector

wg=zeros(length(x),len);

if type==1
    for i=1:len
        funl=@(xi) exp(k*(xi-i+dx));
        funr=@(xi) exp(k*(-xi+i+dx));
        func=@(xi) 1+k*sin((xi-i+dx)*pi/2/dx)/(pi/2/dx);

        xi=x<i-dx;
        xi=xi(xi>0);
        i_min=length(xi);
        xi=x>i+dx;
        xi=xi(xi<1);
        i_max=length(xi);
        wg(1:i_min,i)=funl(x(1:i_min));
        wg(i_min+1:i_max,i)=func(x(i_min+1:i_max));
        wg(i_max+1:end,i)=funr(x(i_max+1:end));
        
        wg(:,i)=wg(:,i)/max(wg(:,i));

        %Just to see the modes
    %     plot(x,wg(i,:).^2);
    %     hold on
    %     plot((i-dx)*[1 1],[0 2],'k');
    %     plot((i+dx)*[1 1],[0 2],'k');
    %     axis tight;
    end
elseif type==0
    for i=1:len
        xi=x<i-.5;
        xi=xi(xi>0);
        i_min=length(xi);
        xi=x>i+.5;
        xi=xi(xi<1);
        i_max=length(xi);
        
        wg(i_min+1:i_max-1,i)=cos((x(i_min+1:i_max-1)-i)*pi).^2;
    end
end

dt_min=min(diff(t));
if dt_min<0
    error('Your timeline doesnt follow causality');
end
t_max=max(t);
t_min=min(t);
n_t=ceil(t_res*(t_max-t_min)/dt_min);
ti=linspace(t_min,t_max,n_t);
surf_mat=zeros(length(x),n_t);

surf_mat_real(:,1)=wg*real(phi(:,1));
surf_mat_imag(:,1)=wg*imag(phi(:,1));
%Create matrix
phi_ind=1;
for i=1:n_t
    if inter==0
        %No interpolation
        if i~=n_t
            %phi_ind=?
            while t(phi_ind+1)<=ti(i)
                phi_ind=phi_ind+1;
            end 
            surf_mat(:,i)=wg*abs(phi(:,phi_ind)).^2;
        else
            surf_mat(:,i)=wg*abs(phi(:,end)).^2;
            surf_mat_real(:,2)=wg*real(phi(:,end));
            surf_mat_imag(:,2)=wg*imag(phi(:,end));
        end
    elseif inter==1
        if i~=n_t
            while t(phi_ind+1)<=ti(i)
                phi_ind=phi_ind+1;
            end 
            %interpolator parameter
            dt=t(phi_ind+1)-t(phi_ind);
            d=1-(t(phi_ind+1)-ti(i))/dt;
            
            gc=phi(:,phi_ind)'*phi(:,phi_ind+1);
            g=acos(gc);
            gc2=phi(:,phi_ind+1)-phi(:,phi_ind)*gc;
            gc2=gc2/sqrt((gc2'*gc2));
            c=(cos(d*g)*phi(:,phi_ind)+sin(d*g)*gc2);
            c=c/sqrt((c'*c));
            surf_mat(:,i)=wg*(conj(c).*c);
        else
            surf_mat(:,i)=wg*abs(phi(:,end)).^2;
            surf_mat_real(:,2)=wg*real(phi(:,end));
            surf_mat_imag(:,2)=wg*imag(phi(:,end));
        end    
    elseif inter==2
        if i~=n_t
            while t(phi_ind+1)<=ti(i)
                phi_ind=phi_ind+1;
            end 
            %interpolator parameter
            dt=t(phi_ind+1)-t(phi_ind);
            d=1-(t(phi_ind+1)-ti(i))/dt;
            
            gc=abs(phi(:,phi_ind))'*abs(phi(:,phi_ind+1));
            g=acos(gc);
            gc2=abs(phi(:,phi_ind+1))-abs(phi(:,phi_ind))*gc;
            gc2=gc2/sqrt((gc2'*gc2));
            c=(cos(d*g)*abs(phi(:,phi_ind))+sin(d*g)*gc2);
            c=c/sqrt((c'*c));
            surf_mat(:,i)=wg*(conj(c).*c);
        else
            surf_mat(:,i)=wg*abs(phi(:,end)).^2;
            surf_mat_real(:,2)=wg*real(phi(:,end));
            surf_mat_imag(:,2)=wg*imag(phi(:,end));
        end
    end
end
end

