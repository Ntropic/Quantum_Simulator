function  gate = gray_gate_substitute( gate,halving,subsi,i )
%adds substitutions phi in gate names for the circuit representation

if halving==0
    if min(abs([0 -1]+abs(mod(sqrt(subsi(i))/2,1)-1)))>10^-14
        if  min(abs([0 -1]+mod(sqrt(subsi(i)),1)))<10^-14
            if abs(round(sqrt(subsi(i)))-1)>10^-14
                gate.circuit_subs{end+1}=[{'\phi',[num2str(round(sqrt(subsi(i)))) '\phi' ]}];
            else
                 gate.circuit_subs{end+1}=[];
            end
        else
            gate.circuit_subs{end+1}=[{'\phi',['\sqrt{' num2str(round(subsi(i))) '}\phi' ]}];
        end
    else
        if round(sqrt(subsi(i))/2)==1
            gate.circuit_subs{end+1}=[{'\frac{\phi}{2}','\phi'}];
        else   
            gate.circuit_subs{end+1}=[{'\frac{\phi}{2}',[num2str(round(sqrt(subsi(i))/2)) '\phi' ]}];
        end
    end
    
else
    subsi2=subsi(i);
    if  min(abs([0 -1]+abs(mod(sqrt(subsi2)/4,1)-1)))>10^-14  %sqrt(subsi2)/4 isn't natural number
        if  min(abs([0 -1]+mod(sqrt(subsi2)/2,1)))<10^-14 || mod(sqrt(subsi2)/2,1)>1-10^-14 %sqrt(subsi2)/2 is natural
            if abs(round(sqrt(subsi2)/2)-1)>10^-14
                gate.circuit_subs{end+1}=[{'\phi',[num2str(round(sqrt(subsi2)/2)) '\phi' ]}];
            else
                gate.circuit_subs{end+1}=[];
            end
        else
            if mod(sqrt(subsi2),1)<10^-14 || mod(sqrt(subsi2),1)>1-10^-14 %subsi2 is natural
                if round(sqrt(subsi2))==1
                    gate.circuit_subs{end+1}=[{'\frac{\phi}{2}',['\frac{\phi}{4}' ]}];
                else
                    gate.circuit_subs{end+1}=[{'\frac{\phi}{2}',['\frac{' num2str(round(sqrt(subsi2))) '}{4}\phi' ]}];
                end
            else
                gate.circuit_subs{end+1}=[{'\frac{\phi}{2}',['\frac{\sqrt{' num2str(round(subsi2)) '}}{4}\phi' ]}];
            end
        end
    else
        gate.circuit_subs{end+1}=[];
    end
end

end

