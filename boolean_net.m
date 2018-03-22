clc
clear all
close all

figure = input('What figure do you wish to generate ? :  ')

switch figure
    case 15
        
        % There are 8 inputs in total
        num_in = 8;
        
        % There are 6 faults considered for this simulation
        num_fault = 6;
        
        % There are 7 drugs considered for this simulation
        num_drugs = 7;
        
        % We consider no inputs, TRAIL and ER Stress input
        
        input = [0,32,128];
        
        % We consider no faults, and a case with DR5 and STAT3 faults
        
        fault = [0,3,3,3];
        
        % We consider no drugs, SH5-07 and CT drug action
        
        drug = [0,0,2,64];
        
        for j = 1 : length(fault)
            
            faults = fault(j);
            drugs = drug(j);
            
            for i = 1 : length(input)
                
                inputs = input(i);
                
                [drug_vector, apoptosis_ratio] = boolean_network_function(num_in,num_fault,num_drugs,inputs,faults,drugs);
                
                apop(i,j) = apoptosis_ratio;
                
            end
        end
        
        % Generating Figure 16
        c = categorical({'00000000','00100000','10000000'});
        xlabel('Inputs')
        ylabel('Apoptosis Ratio')
        bar(c,apop)
        legend('No Faults','DR5 and STAT3 Fault','SH5-07 action on DR5 and STAT3 Fault','CT action on DR5 and STAT3 Fault')
        
    case 16
        
        % There are 8 inputs in total
        num_in = 8;
        
        % There are 6 faults considered for this simulation
        num_fault = 6;
        
        % There are 7 drugs considered for this simulation
        num_drugs = 7;
        
        % The simulation is done with TRAIL as the active input
        % The input vector should be [0010000]
        
        inputs = 32;
        
        % We consider simultaneous presence of all faults
        faults = 63;
        
        % We consider CT in combination with a single drug
        drugs = [64,65,66,68,72,80,96];
        
        [drug_vector, apoptosis_ratio] = boolean_network_function(num_in,num_fault,num_drugs,inputs,faults,drugs);
        
        % Generating Figure 16
        c = categorical({'1000000','1000001','1000010','1000100','1001000','1010000','1100000'});
        xlabel('Cryptotanshinone with one drug at a time')
        ylabel('Net Apoptosis Ratio')
        bar(c,apoptosis_ratio')
        legend('Fault : 111111')
        
    case 17
        
        % There are 8 inputs in total
        num_in = 8;
        
        % There are 6 faults considered for this simulation
        num_fault = 6;
        
        % There are 7 drugs considered for this simulation
        num_drugs = 7;
        
        % The simulation is done with TRAIL as the active input
        % The input vector should be [0010000]
        
        inputs = 32; 
        
        % We consider all possible combination of the faults
                faults = 0 : power(2,num_fault)-1;
        
        % We consider all possible combination of the drugs with
        % Cryptotanshinone
                drugs = power(2,num_drugs-1): power(2,num_drugs)-1;
        
        [drug_vector, apoptosis_ratio] = boolean_network_function(num_in,num_fault,num_drugs,inputs,faults,drugs);
        
                % Generate Results.xlsx and Figure 17
                T = table(drug_vector',apoptosis_ratio');
                filename = 'results_17.xlsx';
                writetable(T,filename,'Sheet',1,'Range','B1')
        
        
    case 18
        
        % There are 8 inputs in total
        num_in = 8;
        
        % There are 6 faults considered for this simulation
        num_fault = 6;
        
        % There are 7 drugs considered for this simulation
        num_drugs = 6;
        
        % The simulation is done with TRAIL as the active input
        % The input vector should be [0010000]
        
        inputs = 32;
        
        % We consider all possible combination of the faults
        faults = 0 : power(2,num_fault)-1;
        
        % We consider all possible combination of the drugs without
        % Cryptotanshinone
        drugs = 0 : power(2,num_drugs)-1;
        
        [drug_vector, apoptosis_ratio] = boolean_network_function_nocryp(num_in,num_fault,num_drugs,inputs,faults,drugs);
        
                % Generate Results.xlsx and Figure 18
                T = table(drug_vector',apoptosis_ratio');
                filename = 'results_18.xlsx';
                writetable(T,filename,'Sheet',1,'Range','B1')
        
        
    otherwise
        disp('Please choose from 15,16,17 or 18.')
end