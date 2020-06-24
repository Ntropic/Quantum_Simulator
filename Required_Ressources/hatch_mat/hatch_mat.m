function [  ] = hatch_mat( ax,mode,mat,density,x,y )
%Create hatchet plot via mat

if length(mode)==0
    mode='any_borders';
end

x=x-0.5;
y=y-0.5;
if length(findstr(mode,'any'))>0
    mat=(mat>0);
    for i=1:size(mat,1)
        for j=1:size(mat,2)
            if mat(i,j)==1
                for k=1:density
                    dx=x(j)+(k-1)/density;
                    dy=y(i)+1-(k-1)/density;
                    plot(ax,[dx x(j)+1],[y(i) dy],'k');
                    plot(ax,[x(j) dx],[dy y(i)+1],'k');
                end
                if length(findstr(mode,'borders'))>0
                    if i>1
                        if mat(i-1,j)==0
                            plot(ax,[x(j) x(j)+1],[y(i) y(i)],'k');
                        end
                    end
                    if i<size(mat,1)
                        if mat(i+1,j)==0
                            plot(ax,[x(j) x(j)+1],[y(i)+1 y(i)+1],'k');
                        end
                    end
                    
                    if j>1
                        if mat(i,j-1)==0
                            plot(ax,[x(j) x(j)],[y(i) y(i)+1],'k');
                        end
                    end
                    if j<size(mat,2)
                        if mat(i,j+1)==0
                            plot(ax,[x(j)+1 x(j)+1],[y(i) y(i)+1],'k');
                        end
                    end
                end
            end
        end
    end
end

end

