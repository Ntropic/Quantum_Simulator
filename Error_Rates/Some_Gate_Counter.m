function [ toffolis, cus ] = Some_Gate_Counter( gate )
%GATE_COUNTER counts single and double qubit gates within the defined gates
%Counts toffolis, CU's and CZ's
gates=gate.steps.gates;

toffolis=0;
cus=0;

for i=1:length(gates)
    if length(findstr(gates{i},'toffoli'))>0
        toffolis=toffolis+1;
    elseif length(findstr(gates{i},'C'))>0
        cus=cus+1;
    end
end
end