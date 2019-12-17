%Script to analyze the performance of the Neural Networks trained with the
%Version 4 of the processed and smoothed experimental data.
%Two different techniques are used to improve generalization: Early Stop
%and Bayesian Regularization
clear all;
close all;
clc;
%This file is required to calculate the R2 value of each network
load('TrainingSetProcStressIndV4.mat'); %PmatData

file_names = {"NN_NatR_Bayesian_Tansig.mat";
    "NN_NatR_Bayesian_Logsig.mat";
    "NN_NatR_Bayesian_Purelin.mat";
    "NN_NatR_LM_Tansig.mat";
    "NN_NatR_LM_Logsig.mat";
    "NN_NatR_LM_Purelin.mat"};

%Obtain dimensions the structure
A = length(file_names); %Number of files
load(file_names{1});
fields.N = fieldnames(neuralnets);
N = length(fields.N); %Number of materials
fields.M = fieldnames(neuralnets.NatR);
M = length(fields.M); %Number of topologies
%Create figures placeholders
%Each figure will have six sublots, one for each training method, showing
%the error against number of neurons
for i=1:M
   fig_performance{i} = figure;
   fig_performance{i}.Name = fields.M{i};
   fig_performance{i}.WindowState = 'maximized';
   fig_nmad{i} = figure;
   fig_nmad{i}.Name = fields.M{i};
   fig_nmad{i}.WindowState = 'maximized';
   fig_r2{i} = figure;
   fig_r2{i}.Name = fields.M{i};
   fig_r2{i}.WindowState = 'maximized';
end

for i=1:A %Number of files
   load(file_names{i});   
   NNs = neuralnets;   
   for j=3 %Number of materials
       for k=1:M %Number of topologies
           figure(fig_performance{k});
           %Performance of network per neuron combination
           T = length( NNs.(fields.N{j}).(fields.M{k}) ); %Total nuerons tested
           perf = 0.*T;
           for L=1:T
            perf(L) = NNs.(fields.N{j}).(fields.M{k})(L).performance.best_perf;
           end           
           subplot(2,3,i);
           semilogy(1:length(perf),perf);
           title(strcat("Method: ", file_names{i} ));
           ylabel("MSE");
           xlabel("Number of Neurons");
           grid on;
           grid minor;           
           %% Calculate and plot the Normalized Mean Absolute Difference (NMAD)
            %and the R2 value
            fields.P = fieldnames(PmatData.(fields.N{j}));
            P = length(fields.P);
            for m=1:P %Number of tests
                nmad = zeros(T,1);
                r2 = zeros(T,1);
                switch fields.P{m} %Type of test
                    case "disR50" 
                        tempRate = 50/60/33; %Specimen length = 33mm
                    case "disR250"
                        tempRate = 250/60/33;
                    case "disR500"
                        tempRate = 500/60/33;
                end
                for n=1:T %Number of neurons
                    b = length(PmatData.(fields.N{j}).(fields.P{m}));
                    rnb = randi(b); %Choose one of the specimens data at random
                    %Experimental Data
                    X_strain = PmatData.(fields.N{j}).(fields.P{m})(rnb).pstrain;
                    X_strain_1 = circshift(X_strain,1,1);
                    X_strain_2 = circshift(X_strain,2,1);
                    X_strainRate = tempRate.*(ones(1,length(X_strain)));
                    X_strainRate_1 = circshift(X_strainRate,1,1);
                    X_strainRate_2 = circshift(X_strainRate,2,1);
                    Y_stress = PmatData.(fields.N{j}).(fields.P{m})(rnb).pstress*1e-6; %Convert to MPa
                    Y_stress_1 = circshift(Y_stress,1,1);
                    Y_stress_2 = circshift(Y_stress,2,1);
                    %Prediction
                    switch k %Network Topologies
                        case 1
                            inputs = X_strain;
                        case 2
                            inputs = [X_strain;X_strain_1];
                        case 3
                            inputs = [X_strain;X_strain_1;Y_stress_1];
                        case 4
                            inputs = [X_strain;X_strainRate];
                        case 5
                            inputs = [X_strain; X_strain_1; X_strainRate];
                        case 6
                            inputs = [X_strain; X_strain_1; X_strainRate;...
                                X_strainRate_1];
                        case 7
                            inputs = [X_strain; X_strain_1; X_strainRate;...
                                X_strainRate_1; Y_stress_1];
                    end
                    y_stress = NNs.(fields.N{j}).(fields.M{k})(n).NN( inputs );
                    %NMAD
                    nmad(n) = 100*sum ( abs(y_stress-Y_stress) ) ./ sum( abs(Y_stress) );
                    %Coefficient of Determination
                    r2(n) = 1 - sum( (y_stress-Y_stress).^2 )./ sum( ( Y_stress - mean(Y_stress)).^2 );
                end
                %Plots for NMAD
                figure(fig_nmad{k})
                hold on
                subplot(2,3,i);
                plot(nmad)
                title(strcat("Method: ", file_names{i} ));
                ylabel("Normalized Mean Absolute Difference");
                xlabel("Number of Neurons");
                ylim([0 100]);
                grid on;
                grid minor;
                hold off
                
                %Plots for R^2
                figure(fig_r2{k})
                hold on
                subplot(2,3,i);
                plot(r2)
                title(strcat("Method: ", file_names{i} ));
                ylabel("Coefficient of Determination (R^2)");
                xlabel("Number of Neurons");
                ylim([0 1]);
                grid on;
                grid minor;
                hold off
            end
       end
   end  
end


           
