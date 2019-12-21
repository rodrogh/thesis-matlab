%Script to analyze the performance of the Neural Networks trained with the
%Version 4 of the processed and smoothed experimental data.
%Two different techniques are used to improve generalization: Early Stop
%and Bayesian Regularization
clear all;
close all;
clc;
%Append working directory due to GIT integration
dir = "../Experimental Work/";
%This file is required to calculate the R2 value of each network
load( strcat(dir,"TrainingSetProcStressIndV4.mat")); %PmatData
load(strcat(dir,"FormattedTrainingSetIndV4.mat")); %trainData & validationData


file_names = [ strcat(dir,"NN_NatR_Bayesian_Tansig.mat");
    strcat(dir,"NN_NatR_Bayesian_Logsig.mat");
    strcat(dir,"NN_NatR_Bayesian_Purelin.mat");
    strcat(dir,"NN_NatR_LM_Tansig.mat");
    strcat(dir,"NN_NatR_LM_Logsig.mat");
    strcat(dir,"NN_NatR_LM_Purelin.mat")];

%Obtain dimensions the structure
A = length(file_names); %Number of files
load( strcat("../Experimental Work/" ,file_names(1)));
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
%     fig_mse_u{i} = figure;
%     fig_mse_u{i}.Name = fields.M{i};
%     fig_mse_u{i}.WindowState = 'maximized';
%     fig_r2{i} = figure;
%     fig_r2{i}.Name = fields.M{i};
%     fig_r2{i}.WindowState = 'maximized';
end

for i=1:A %Number of files
    load(file_names(i));
    NNs = neuralnets;
    for j=3 %Number of materials - Does not work for rubber bands
        %Training Set (Transpose to match NN() requirements)
        X_strain = trainData.(fields.N{j}).input_stress(:,1)';
        X_strain_1 = trainData.(fields.N{j}).input_stress(:,2)';
        X_strain_2 = trainData.(fields.N{j}).input_stress(:,3)';
        X_strainRate = trainData.(fields.N{j}).input_stress(:,4)';
        X_strainRate_1 = trainData.(fields.N{j}).input_stress(:,5)';
        X_strainRate_2 = trainData.(fields.N{j}).input_stress(:,6)';
        Y_stress = trainData.(fields.N{j}).output_stress(:,1)'.*1e-6; %Convert to MPa
        Y_stress_1 = trainData.(fields.N{j}).output_stress(:,2)'.*1e-6;
        Y_stress_2 = trainData.(fields.N{j}).output_stress(:,3)'.*1e-6;
        %Unknown Set
        Xu_strain = validationData.(fields.N{j}).input_stress(:,1)';
        Xu_strain_1 = validationData.(fields.N{j}).input_stress(:,2)';
        Xu_strain_2 = validationData.(fields.N{j}).input_stress(:,3)';
        Xu_strainRate = validationData.(fields.N{j}).input_stress(:,4)';
        Xu_strainRate_1 = validationData.(fields.N{j}).input_stress(:,5)';
        Xu_strainRate_2 = validationData.(fields.N{j}).input_stress(:,6)';
        Yu_stress = validationData.(fields.N{j}).output_stress(:,1)'.*1e-6; %Convert to MPa
        Yu_stress_1 = validationData.(fields.N{j}).output_stress(:,2)'.*1e-6;
        Yu_stress_2 = validationData.(fields.N{j}).output_stress(:,3)'.*1e-6;
        for k=1:M %Number of topologies
            %% Local Variables Initialization
            T = length( NNs.(fields.N{j}).(fields.M{k}) ); %Total nuerons tested
            neuronsNr = 1:T; %known data
            best_perf = 0.*T;
            u_perf = 0.*T; %unknown data
            train_R2 = 0.*T; %R^2
            val_R2 = 0.*T;
            test_R2 = 0.*T;
            %Selection of network topology 
            switch k %Network Topologies
                case 1
                    inputs = X_strain;
                    inputs_u = Xu_strain;
                case 2
                    inputs = [X_strain;X_strain_1];
                    inputs_u = [Xu_strain;Xu_strain_1];
                case 3
                    inputs = [X_strain;X_strain_1;Y_stress_1];
                    inputs_u = [Xu_strain;Xu_strain_1;Yu_stress_1];
                case 4
                    inputs = [X_strain;X_strainRate];
                    inputs_u = [Xu_strain;Xu_strainRate];
                case 5
                    inputs = [X_strain; X_strain_1; X_strainRate];
                    inputs_u = [Xu_strain; Xu_strain_1; Xu_strainRate];
                case 6
                    inputs = [X_strain; X_strain_1; X_strainRate;...
                        X_strainRate_1];
                    inputs_u = [Xu_strain; Xu_strain_1; Xu_strainRate;...
                        Xu_strainRate_1];
                case 7
                    inputs = [X_strain; X_strain_1; X_strainRate;...
                        X_strainRate_1; Y_stress_1];
                    inputs_u = [Xu_strain; Xu_strain_1; Xu_strainRate;...
                        Xu_strainRate_1; Yu_stress_1];
            end
            %% Performance and r^2 Extraction
            for n=1:T
                %% MSE for Known data
                best_perf(n) = NNs.(fields.N{j}).(fields.M{k})(n).performance.best_perf;
                %% MSE of Unknown data
                y_stress_u = NNs.(fields.N{j}).(fields.M{k})(n).NN( inputs_u );
                u_perf(n) = mean(y_stress_u - Yu_stress ).^2;
                %% R2 Coefficient during training, validation, testing
                %Use indices found in performance struct in combination with the trainData variable
                trainInd = NNs.(fields.N{j}).(fields.M{k})(n).performance.trainInd;
                valInd = NNs.(fields.N{j}).(fields.M{k})(n).performance.valInd;
                testInd = NNs.(fields.N{j}).(fields.M{k})(n).performance.testInd;
                %Network Prediction
                y_stress_train = NNs.(fields.N{j}).(fields.M{k})(n).NN( inputs(:,trainInd) );
                y_stress_val = NNs.(fields.N{j}).(fields.M{k})(n).NN( inputs(:,valInd) );
                y_stress_test = NNs.(fields.N{j}).(fields.M{k})(n).NN( inputs(:,testInd) );
                %Coefficient of Determination
                train_R2(n) = 1 - sum( (y_stress_train-Y_stress(trainInd)).^2 )./ sum( ( Y_stress(trainInd) - mean(Y_stress(trainInd))).^2 );
                val_R2(n) = 1 - sum( (y_stress_val-Y_stress(valInd)).^2 )./ sum( ( Y_stress(valInd) - mean(Y_stress(valInd))).^2 );
                test_R2(n) = 1 - sum( (y_stress_test-Y_stress(testInd)).^2 )./ sum( ( Y_stress(testInd) - mean(Y_stress(testInd))).^2 );
            end
            
            %% Plots
            %MSE Left
            figure(fig_performance{k});
            subplot(2,3,i);
            yyaxis left
            semilogy(neuronsNr,best_perf,neuronsNr,u_perf);
            ylabel("Mean Squared Error (MSE)");
            ylim([10e-6 10]);
            %R2 Right
            yyaxis right
            if contains( file_names{i}, 'Bayesian')
                plot(neuronsNr,train_R2,neuronsNr,test_R2);
                legend('Training Dataset','Unknown Dataset',...%Left
                'Training','Testing','Location','best'); %Right
            else
                plot(neuronsNr,train_R2,neuronsNr,val_R2,neuronsNr,test_R2);
                legend('Training Dataset','Unknown Dataset',...%Left
                'Training','Validation','Testing','Location','best'); %Right
            end
            ylabel("Coefficient of Determination (R^2)");
            ylim([0 1]);
            title(strcat("Method: ", file_names(i) ));
            xlabel("Number of Neurons");
            grid on;
            grid minor;
            
%             %Plots for R^2
%             figure(fig_r2{k})
%             subplot(2,3,i);
%             title(strcat("Method: ", file_names(i) ));
%             legend();
%             xlabel("Number of Neurons");
%             ylim([0 1]);
%             grid on;
%             grid minor;                      
        end
    end
end



