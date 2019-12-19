%Script to observe the impact of increasing the number of nuerons on a
%feedforward ANN. The matlab tables with the formatted datsets must be
%loaded
function ANN_Neuron_Test(fname,trainFcn,activationFcn,functionality)
% clear all;
close all;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Initialization
load('../Experimental work/FormattedTrainingSetInd.mat');
%Defining number of materials
N = length( fieldnames(trainData) ); %Value of N is valid for all cases
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(trainData);
T = 25; %Maximum number of neurons to test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Behaviour variables
saveData = true;
%Choose between: 'single training' 'multi training' 'getPrediction' 'plots'
%'getError'
%If Analysis is desired, shoose input type: 'unknown' 'known'
inputType = 'known';
%Neural network predictions take too much time to compute, if only
%the prediction of the mathematical model is desired, choose
analysisType = 'complete';

switch functionality
    %%
    case 'single training'
        disp( 'Test Stage: START of Single Layer Analysis' );
        neuralnets = trainData;
        for i=1:N
            M = length( fieldnames(trainData.(fields.N{i})) );
            %Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
            fields.M = fieldnames(trainData.(fields.N{i}));
            
            for j=1:M
                neuralnets.(fields.N{i}).(fields.M{j}) = [];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %First configuration - FeedForward Single Layer
                %Input: Displacement
                %Output: Load
                x = trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';
                t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                [NN, NNPerform] = ffnn_single(x,t,T);
                neuralnets.(fields.N{i}).(fields.M{j}).ff1 = struct('NN',NN,'performance',NNPerform);
                
                %Second configuration - FeedForward Single Layer
                %Inputs: Displacement - Displacement(t-1)
                %Output: Load
                x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)'];
                t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                [NN, NNPerform] = ffnn_single(x,t,T);
                neuralnets.(fields.N{i}).(fields.M{j}).ff2 = struct('NN',NN,'performance',NNPerform);
                
                %Third configuration - FeedForward Single Layer
                %Inputs: Displacement - Displacement(t-1) - Displacement(t-2)
                %Output: Load
                x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,3)'];
                t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                [NN, NNPerform] = ffnn_single(x,t,T);
                neuralnets.(fields.N{i}).(fields.M{j}).ff3 = struct('NN',NN,'performance',NNPerform);
                
                %Fourth configuration - FeedForward Single Layer
                %Inputs: Displacement - Displacement(t-1) -
                %Displacement(t-2) - Load(t-1)
                %Output: Load
                x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,3)';...
                    trainData.(fields.N{i}).(fields.M{j}).outputs(:,2)'];
                t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                [NN, NNPerform] = ffnn_single(x,t,T);
                neuralnets.(fields.N{i}).(fields.M{j}).ff4 = struct('NN',NN,'performance',NNPerform);
            end
            disp( strcat( 'Test Stage: ', fields.N{i} ) );
        end
        disp( 'Test Stage: FINISH of Single Layer Analysis' );
        %After analysis of a single layer with T nuerons, Extraction of the
        %NN with best performance
        nn_best_single = getBestNN(neuralnets);
        
        disp('Done!');
        if saveData
            disp('Saving Data Single...')
            save('NN_Best.mat','nn_best_single');
        end
        %%
    case 'single trainingV2'
        disp( 'Test Stage: START of Single Layer Analysis V4' );
        %Initialization
        load('../Experimental work/FormattedTrainingSetIndV4.mat');
        %Defining number of materials
        N = length( fieldnames(trainData) ); %Value of N is valid for all cases
        %Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
        fields.N = fieldnames(trainData);
        T = 30; %Maximum number of neurons to test
        neuralnets = trainData;
        for i=3
            %Assigning inputs and outputs to new variables
            inputs = trainData.(fields.N{i}).input_stress;
            outputs = trainData.(fields.N{i}).output_stress*1e-6; %Convert to MPa
            disp( strcat( 'Test Stage: ', fields.N{i} ) );
            neuralnets.(fields.N{i}) = [];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Rate Independent Neural Network Acrchitectures
            %First configuration - FeedForward Single Layer
            %Input: Displacement
            %Output: Load
            x = inputs(:,1)';
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff1 = struct('NN',NN,'performance',NNPerform);
            
            %Second configuration - FeedForward Single Layer
            %Inputs: Displacement - Displacement(t-1)
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,2)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff2 = struct('NN',NN,'performance',NNPerform);
            
            %Third configuration - FeedForward Single Layer
            %Inputs: Displacement - Displacement(t-1) - Load(t-1)
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,2)';...
                outputs(:,2)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff3 = struct('NN',NN,'performance',NNPerform);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% Rate Dependent Neural Network Architectures
            %Fourth configuration - FeedForward Single Layer
            %Inputs: Displacement - DisRate
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,4)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff4 = struct('NN',NN,'performance',NNPerform);
            
            %Fifth configuration - FeedForward Single Layer
            %Inputs: Displacement - Displacement(t-1) - DisRate
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,2)';...
                inputs(:,4)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff5 = struct('NN',NN,'performance',NNPerform);
            
            %Sixth configuration - FeedForward Single Layer
            %Inputs: Displacement - Displacement(t-1) - DisRate -
            %DisRate(t-1)
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,2)';...
                inputs(:,4)';...
                inputs(:,5)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff6 = struct('NN',NN,'performance',NNPerform);
            
            %Seventh configuration - FeedForward Single Layer
            %Inputs: Displacement - Displacement(t-1) - DisRate -
            %DisRate(t-1) - Load(t-1)
            %Output: Load
            x = [ inputs(:,1)';...
                inputs(:,2)';...
                inputs(:,4)';...
                inputs(:,5)';...
                outputs(:,2)'];
            t = outputs(:,1)';
            [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn);
            neuralnets.(fields.N{i}).ff7 = struct('NN',NN,'performance',NNPerform);
        end
        disp( 'Test Stage: FINISH of Single Layer Analysis' );
        %After analysis of a single layer with T neurons, Extraction of the
        %NN with best performance
        %         load('NN_BestV2NatR_RI.mat');
        %         nn_best_single = getBestNN(nn_best_single);
        %getBestNN is commented for Version 3 because I need to look
        %carefully at the performance of all the trained networks
        %         nn_best_single = getBestNN(neuralnets);
        disp('Done!');
        if saveData
            disp('Saving Data Single...')
            %             save('NN_BestV3NatR.mat','nn_best_single');
            fileName = strcat('../Experimental work',fname,'.mat');
            save(fileName,'neuralnets');
        end
    case 'multi training'
        disp( 'Test Stage: START of Multi Layer Analysis' );
        load('NN_BestV3NatR.mat'); %nn_best_single
        for i=3
            %Test variations (disR)
            M = length( fieldnames(trainData.(fields.N{i})) );
            fields.M = fieldnames(trainData.(fields.N{i}));
            
            for j=1:M
                %                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %                 %First config variant - FeedForward Two Layers
                %                 %Input: Displacement
                %                 %Output: Load
                %                 x = trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';
                %                 t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                %                 best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff1.index;
                %
                %                 [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                %                 nn_best_single.(fields.N{i}).(fields.M{j}).ff1_m = struct('NN',NN,'performance',NNPerform);
                %
                %                 %Second config variant - FeedForward Two Layers
                %                 %Inputs: Displacement - Displacement(t-1)
                %                 %Output: Load
                %                 x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)'];
                %                 t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                %                 best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff2.index;
                %
                %                 [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                %                 nn_best_single.(fields.N{i}).(fields.M{j}).ff2_m = struct('NN',NN,'performance',NNPerform);
                %
                %                 %Third config - FeedForward Two Layers
                %                 %Inputs: Displacement - Displacement(t-1) - Displacement(t-2)
                %                 %Output: Load
                %                 x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).inputs(:,3)'];
                %                 t = trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                %                 best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff3.index;
                %
                %                 [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                %                 nn_best_single.(fields.N{i}).(fields.M{j}).ff3_m = struct('NN',NN,'performance',NNPerform);
                %
                %                 %Fourth config - FeedForward Two Layers
                %                 %Inputs: Displacement - Displacement(t-1) - ...
                %                 %Displacement(t-2) - Load(t-1)
                %                 %Output: Load
                %                 x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).inputs(:,3)';...
                %                     trainData.(fields.N{i}).(fields.M{j}).outputs(:,2)'];
                %                 t =  trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                %                 best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff4.index;
                %
                %                 [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                %                 nn_best_single.(fields.N{i}).(fields.M{j}).ff4_m = struct('NN',NN,'performance',NNPerform);
                %
                %Sixth config - FeedForward Two Layers
                %Inputs: Displacement - Displacement(t-1)
                %Output: Load
                x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)'];
                t =  trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff6.index;
                
                [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                nn_best_single.(fields.N{i}).(fields.M{j}).ff6_m = struct('NN',NN,'performance',NNPerform);
                
                %Seventh config - FeedForward Two Layers
                %Inputs: Displacement - Displacement(t-1) - Displacement(t-2) ...
                %Output: Load
                x = [ trainData.(fields.N{i}).(fields.M{j}).inputs(:,1)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,2)';...
                    trainData.(fields.N{i}).(fields.M{j}).inputs(:,3)'];
                t =  trainData.(fields.N{i}).(fields.M{j}).outputs(:,1)';
                best_i = nn_best_single.(fields.N{i}).(fields.M{j}).ff7.index;
                
                [NN, NNPerform] = ffnn_single(x,t,T,best_i);
                nn_best_single.(fields.N{i}).(fields.M{j}).ff7_m = struct('NN',NN,'performance',NNPerform);
            end
            disp( strcat( 'Test Stage: ', fields.N{i} ) );
        end
        disp( 'Test Stage: FINISH of Multi Layer Analysis' );
        %This structure contains the data from nn_best_single as well, plus the
        %best neural networks with two layers
        nn_best_multi = getBestNN(nn_best_single);
        
        disp('Done!');
        if saveData
            disp('Saving Data Multi...')
            save('NN_Multi_BestV2_NatR_RI.mat','nn_best_multi');
        end
        %%
    case 'multi trainingV2'
        disp( 'Test Stage: START of Multi Layer Analysis' );
        load('NN_BestV3NatR.mat'); %nn_best_single
        load('FormattedTrainingSetIndV3.mat');
        %Defining number of materials
        N = length( fieldnames(trainData) ); %Value of N is valid for all cases
        %Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
        fields.N = fieldnames(trainData);
        T = 25; %Maximum number of neurons to test
        for i=3
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %First config variant - FeedForward Two Layers
            %Input: Displacement
            %Output: Load
            x = [ trainData.(fields.N{i}).inputs(:,1)';...
                trainData.(fields.N{i}).inputs(:,4)'];
            t = trainData.(fields.N{i}).outputs(:,1)';
            best_i = nn_best_single.(fields.N{i}).ff1.index;
            
            
            
            
            %             %Sixth config - FeedForward Two Layers
            %             %Inputs: Displacement - Displacement(t-1)
            %             %Output: Load
            %             x = [ trainData.(fields.N{i}).inputs(:,1)';...
            %                 trainData.(fields.N{i}).inputs(:,2)'];
            %             t =  trainData.(fields.N{i}).outputs(:,1)';
            %             best_i = nn_best_single.(fields.N{i}).ff6.index;
            %
            %             [NN, NNPerform] = ffnn_single(x,t,T,best_i);
            %             nn_best_single.(fields.N{i}).ff6_m = struct('NN',NN,'performance',NNPerform);
            %
            %             %Seventh config - FeedForward Two Layers
            %             %Inputs: Displacement - Displacement(t-1) - Displacement(t-2) ...
            %             %Output: Load
            %             x = [ trainData.(fields.N{i}).inputs(:,1)';...
            %                 trainData.(fields.N{i}).inputs(:,2)';...
            %                 trainData.(fields.N{i}).inputs(:,3)'];
            %             t =  trainData.(fields.N{i}).outputs(:,1)';
            %             best_i = nn_best_single.(fields.N{i}).ff7.index;
            %
            %             [NN, NNPerform] = ffnn_single(x,t,T,best_i);
            %             nn_best_single.(fields.N{i}).ff7_m = struct('NN',NN,'performance',NNPerform);
            
            disp( strcat( 'Test Stage: ', fields.N{i} ) );
        end
        disp( 'Test Stage: FINISH of Multi Layer Analysis' );
        %This structure contains the data from nn_best_single as well, plus the
        %best neural networks with two layers
        %         load('miderror.mat');
        nn_best_multi = getBestNN(nn_best_single);
        
        disp('Done!');
        if saveData
            disp('Saving Data Multi...')
            save('NN_Multi_BestV2_NatR_RD.mat','nn_best_multi');
        end
        %%
    case 'getPrediction'
        tic
        %Initial data
        load('NN_Best.mat'); %nn_best_single
        load('NN_Multi_Best.mat'); %nn_best_multi
        load('PiecewiseLinearization.mat'); %modelFit
        load('TrainingSetProc.mat'); %Raw data
        
        %Placeholder for prediction of NN and MathMod
        predictions = nn_best_multi;
        %Load prediction data to overwrite only mathmod part. Comment NN
        %for loops
        %         load('Prediction_KnownDisR.mat');
        for i=1:N
            %Test variations (disR)
            M = length( fieldnames(nn_best_multi.(fields.N{i})) );
            fields.M = fieldnames(nn_best_multi.(fields.N{i}));
            %% Prediction analysis
            for j=1:M
                P = length( fieldnames(nn_best_multi.(fields.N{i}).(fields.M{j}) ) );
                fields.P = fieldnames(nn_best_multi.(fields.N{i}).(fields.M{j}) );
                %% Input Type
                %Define input signal to test (displacement in milimeters)
                sampleSize = 100;
                switch(inputType)
                    case 'known'
                        file_name = strcat( 'Prediction_KnownDisR_Short.mat');
                        %Testing bit, maybe switch case is unnecesary for
                        %known data
                        tempDis = PmatData.(fields.N{i}).(fields.M{j}).rdis;
                        tempL = length(tempDis);
                        tempIndex = nearest( linspace(1,tempL,sampleSize) );
                        dis = tempDis( tempIndex );
                        time = PmatData.(fields.N{i}).(fields.M{j}).rtime( tempIndex );
                        dt = time(2) - time(1);
                    case 'unknown'
                        disR = 500; %mm/min
                        [dis,dt] = defineInput(disR,0,400,sampleSize);
                        file_name = strcat( 'Prediction_UnknownDisR_', int2str(disR) ,'.mat');
                end
                %% Neural Networks Prediction
                %Obtain prediction of neural networks (single and
                %multi)
                for k=1:P
                    predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) = [];
                    predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) = ...
                        nn_eval( dis, (fields.P{k}), ...
                        nn_best_multi.(fields.N{i}).(fields.M{j}).(fields.P{k}).NN );
                end
                %% Math model prediction
                %Main variable is modelFit
                %Ergonomic Variables
                PL_ke = modelFit.(fields.N{i}).(fields.M{j}).PL_ke;
                disSeg = modelFit.(fields.N{i}).(fields.M{j}).disSeg;
                km = modelFit.(fields.N{i}).(fields.M{j}).km;
                thao = modelFit.(fields.N{i}).(fields.M{j}).thao;
                %Starts from 2 because of initial values set to zero
                force = 0.*dis;
                %The method use past values. Therefore the very first value
                %must be initialized
                force(1) = (( PL_ke(1) + km) * dis(1)  + ...
                    PL_ke(1)*dis(1)*dt/thao ) / (1 + dt/thao);
                %Calculation of remainder elements
                for m=2:length(dis)
                    temp = find( disSeg <= dis(m));
                    if length(temp) > 1
                        k_index = temp(end-1);
                    else
                        k_index = 1;
                    end
                    %Piecewise Linearization Equation for SLS model
                    force(m) = (( sum(PL_ke(1:k_index)) + km) * (dis(m) - dis(m-1)) + ...
                        sum(PL_ke(1:k_index))*dis(m)*dt/thao + force(m-1)) / (1 + dt/thao);
                end
                predictions.(fields.N{i}).(fields.M{j}).model = force;
            end
        end
        
        save(file_name,'predictions','dis');
        toc
        %%
    case 'getPredictionV2'
        %Initial data
        load('NN_BestV2NatRubber.mat'); %nn_best_single
        %         load('NN_Multi_Best.mat'); %nn_best_multi
        load('PiecewiseLinearization.mat'); %modelFit
        load('TrainingSetProc.mat'); %Raw data
        
        %Placeholder for prediction of NN and MathMod
        predictions = nn_best_single;
        %         predictions = nn_best_multi;
        %Load prediction data to overwrite only mathmod part. Comment NN
        %for loops
        %         load('Prediction_KnownDisR.mat');
        for i=3
            %Test variations (disR)
            M = length( fieldnames(PmatData.(fields.N{i})) );
            fields.M = fieldnames(PmatData.(fields.N{i}));
            %% Prediction analysis
            predictions.(fields.N{i}) = [];
            for j=1:M
                %% Input Type
                %Define input signal to test (displacement in milimeters)
                sampleSize = 100;
                switch(inputType)
                    case 'known'
                        file_name = strcat( 'Prediction_KnownDisR_NatR_Short.mat');
                        %Testing bit, maybe switch case is unnecesary for
                        %known data
                        tempDis = PmatData.(fields.N{i}).(fields.M{j}).rdis;
                        tempL = length(tempDis);
                        tempIndex = nearest( linspace(1,tempL,sampleSize) );
                        dis = tempDis( tempIndex );
                        time = PmatData.(fields.N{i}).(fields.M{j}).rtime( tempIndex );
                        dt = time(2) - time(1);
                        disR = time.*0;
                        disR(2:end) = diff(dis)./diff(time);
                    case 'unknown'
                        %This might not work here
                        disR = 500; %mm/min
                        [dis,dt] = defineInput(disR,0,400,sampleSize);
                        file_name = strcat( 'Prediction_UnknownDisR_NatR_', int2str(disR) ,'.mat');
                end
                %% Neural Networks Prediction
                %Obtain prediction of neural networks
                P = length( fieldnames(nn_best_single.(fields.N{i})) );
                fields.P = fieldnames(nn_best_single.(fields.N{i}));
                
                for k=1:P
                    predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) = ...
                        nn_evalV2( dis, disR, fields.P{k}, ...
                        nn_best_single.(fields.N{i}).(fields.P{k}).NN );
                end
            end
        end
        if saveData
            save(file_name,'predictions');
        end
        %%
    case 'getError'
        pfile = 'KnownDisR.mat';
        load( strcat('Prediction_',pfile) );
        %Raw data
        load('TrainingSetProc.mat');
        %Create placeholder
        errors = predictions;
        for i=1:N
            %Test variations (disR)
            M = length( fieldnames(predictions.(fields.N{i})) );
            fields.M = fieldnames(predictions.(fields.N{i}));
            %% Error Extraction
            for j=1:M
                
                figure('Name',strcat(fields.N{i},' - ',fields.M{j}),'NumberTitle','off');
                
                P = length( fieldnames(predictions.(fields.N{i}).(fields.M{j}) ) );
                fields.P = fieldnames(predictions.(fields.N{i}).(fields.M{j}) );
                plotData = zeros(P,5);
                for k=1:P
                    %Prediction - Target (Error)
                    pt = predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) - ...
                        PmatData.(fields.N{i}).(fields.M{j}).pload;
                    %Mean Target  - Target (Variation)
                    mt = mean( PmatData.(fields.N{i}).(fields.M{j}).pload ) - ...
                        PmatData.(fields.N{i}).(fields.M{j}).pload ;
                    %Root Mean Squared Error
                    rmse(k) = sqrt( mean(pt.^2) );
                    %Relative Squared Error
                    rse(k) = sum(pt.^2)/sum(mt.^2);
                    %Mean Absolute Error
                    mae = mean( abs(pt) );
                    %Relative Absolute Error
                    rae = sum( abs(pt) )/sum( abs(mt) );
                    %Coefficient of Determination (R-Squared)
                    r2 = 1 - ( sum( pt.^2) / sum( mt.^2));
                    %Saving Data into the Structure
                    errors.(fields.N{i}).(fields.M{j}).(fields.P{k}) = ...
                        struct('rmse',rmse,'rse',rse,'mae',mae,'rae',rae,'r2',r2);
                    %                     errors(k,:) = [ rmse, rse, mae, rae, r2];
                    
                    %Plot for all 5 error methods
                    %                     subplot(3,3,k);
                    %                     bar( [rmse rse mae rae r2] );
                    %                     title( fields.P{k} );
                    %                     xticks( 1:5 );
                    %                     xticklabels( {'rmse';'rse';'mae';'rae';'r2' } );
                end
                
                bar(rse);
                xticks( 1:P );
                xticklabels( fields.P );
            end
        end
        %Saving into data file
        save( strcat('Errors_',pfile),'errors' );
        
    case 'plots'
        load('PiecewiseLinearization.mat'); %For plotting purpouses
        switch(inputType)
            case 'known'
                load('Prediction_KnownDisR.mat');
            case 'unknown'
                load('Prediction_UnknownDisR.mat');
        end
        for i=1:N
            %Test variations (disR)
            M = length( fieldnames(predictions.(fields.N{i})) );
            fields.M = fieldnames(predictions.(fields.N{i}));
            for j=1:M
                %Topologies variation (ff1...)
                P = length( fieldnames(predictions.(fields.N{i}).(fields.M{j}) ) );
                fields.P = fieldnames(predictions.(fields.N{i}).(fields.M{j}) );
                
                figure('Name',strcat(fields.N{i},' - ',fields.M{j}),'NumberTitle','off');
                
                ax1 = subplot(1,2,1);
                ax2 = subplot(1,2,2);
                hold(ax1,'on')
                hold(ax2,'on')
                for k=1:P
                    %                     subplot((P-1)/2,2,k);
                    if endsWith( fields.P{k}, '_m')
                        plot(ax2,dis,predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}),':',...
                            'DisplayName',(fields.P{k}));
                    else
                        if strcmp(fields.P{k},'model')
                            plot(ax1,dis,predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}),...
                                'DisplayName',(fields.P{k}));
                            plot(ax2,dis,predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}),...
                                'DisplayName',(fields.P{k}));
                        else
                            plot(ax1,dis,predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}),':',...
                                'DisplayName',(fields.P{k}));
                        end
                    end
                end
                xlim(ax1,[0 modelFit.(fields.N{i}).(fields.M{j}).uDis]);
                xlim(ax2,[0 modelFit.(fields.N{i}).(fields.M{j}).uDis]);
                xlabel(ax1,'Dis');
                ylabel(ax1,'Load');
                legend(ax1);
                legend(ax2);
                hold(ax1,'off')
                hold(ax2,'off')
            end
        end
    case 'plotsV2'
        load('TrainingSetProc.mat');
        known = load('Prediction_KnownDisR_NatR_Short.mat');
        for i=3 %Material
            M = length( fieldnames(known.predictions.(fields.N{i})) );
            fields.M = fieldnames(known.predictions.(fields.N{i}));
            f = figure('Name',strcat(fields.N{i},' New Comparison Analysis'),'NumberTitle','off');
            for j=1:M  %Type of test
                P = length( fieldnames(known.predictions.(fields.N{i}).(fields.M{j} ) ) );
                fields.P = fieldnames(known.predictions.(fields.N{i}).(fields.M{j} ));
                for k=1:P %Type of NN
                    disp(j + (k-1));
                    ax1 = subplot(M,P, k + P*(j-1));
                    hold(ax1,'on');
                    sampleSize = length( known.predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) );
                    tempDis = PmatData.(fields.N{i}).(fields.M{j}).rdis;
                    tempL = length(tempDis);
                    tempIndex = nearest( linspace(1,tempL,sampleSize) );
                    dis = tempDis( tempIndex );
                    plot(ax1, PmatData.(fields.N{i}).(fields.M{j}).rdis (tempIndex ) , ...
                        PmatData.(fields.N{i}).(fields.M{j}).pload (tempIndex ), ...
                        'DisplayName', 'ExpData');
                    plot(ax1, PmatData.(fields.N{i}).(fields.M{j}).rdis (tempIndex ) , ...
                        known.predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}), ...
                        'DisplayName', strcat( 'NN ', upper( fields.P{k} ) ));
                    hold(ax1,'off');
                end
            end
        end
        legend(ax1);
    case 'plots_esp'
        load('PiecewiseLinearization.mat'); %For plotting purpouses
        load('TrainingSetProc.mat');
        known = load('Prediction_KnownDisR_Short.mat');
        u50 = load('Prediction_UnknownDisR_50.mat');
        u150 = load('Prediction_UnknownDisR_150.mat');
        u250 = load('Prediction_UnknownDisR_250.mat');
        u400 = load('Prediction_UnknownDisR_400.mat');
        u500 = load('Prediction_UnknownDisR_500.mat');
        for i=1:N
            %Test variations (disR)
            M = length( fieldnames(known.predictions.(fields.N{i})) );
            fields.M = fieldnames(known.predictions.(fields.N{i}));
            f = figure('Name',strcat(fields.N{i},' Comparison Analysis'),'NumberTitle','off');
            for j=1:M
                %Topologies variation (ff1...)
                P = length( fieldnames(known.predictions.(fields.N{i}).(fields.M{j}) ) );
                fields.P = fieldnames(known.predictions.(fields.N{i}).(fields.M{j}) );
                %Get Best Error per escenario: Known, Unknown, Cross
                %% Known Case
                [e,P_index,sample_i] = getBestError('known',known.predictions,PmatData,i,j,fields);
                ax1 = subplot(M,3, 1 + 3*(j-1));
                hold(ax1,'on')
                plot(ax1, PmatData.(fields.N{i}).(fields.M{j}).rdis (sample_i ) , ...
                    PmatData.(fields.N{i}).(fields.M{j}).pload (sample_i ), ...
                    'DisplayName', 'ExpData');
                plot(ax1, PmatData.(fields.N{i}).(fields.M{j}).rdis(sample_i ) , ...
                    known.predictions.(fields.N{i}).(fields.M{j}).(fields.P{P_index}), ...
                    'DisplayName', strcat( 'BestNN-', upper( fields.P{P_index} ), ' @', num2str(e(P_index),2) ));
                plot(ax1, PmatData.(fields.N{i}).(fields.M{j}).rdis(sample_i ) , ...
                    known.predictions.(fields.N{i}).(fields.M{j}).(fields.P{end}), ...
                    'DisplayName', strcat('Model @',num2str(e(end),2 )));
                xlim(ax1,[0 modelFit.(fields.N{i}).(fields.M{j}).uDis]);
                xlabel(ax1,'Dis');
                ylabel(ax1,'Load');
                legend(ax1);
                title(ax1,strcat( fields.M{j}, 'Known Input') );
                hold(ax1,'off')
                %% Unkown Case
                [e,P_index,~] = getBestError('unknown',u150.predictions,PmatData,i,j,fields);
                ax2 = subplot(M,3, 2 + 3*(j-1));
                hold(ax2,'on')
                plot(ax2, u150.dis , ...
                    u150.predictions.(fields.N{i}).(fields.M{j}).model, ...
                    'DisplayName', 'Model 150');
                plot(ax2, u150.dis , ...
                    u150.predictions.(fields.N{i}).(fields.M{j}). ( fields.P{P_index} ),'--', ...
                    'DisplayName', strcat( 'BestNN 150-', upper( fields.P{P_index} ), ' @', num2str(e(P_index),2) ));
                
                [e,P_index,~] = getBestError('unknown',u400.predictions,PmatData,i,j,fields);
                plot(ax2, u400.dis , ...
                    u400.predictions.(fields.N{i}).(fields.M{j}).model, ...
                    'DisplayName', 'Model 400');
                plot(ax2, u400.dis , ...
                    u400.predictions.(fields.N{i}).(fields.M{j}). ( fields.P{P_index} ),'--', ...
                    'DisplayName', strcat( 'BestNN 400-', upper( fields.P{P_index} ), ' @', num2str(e(P_index),2) ));
                xlim(ax2,[0 modelFit.(fields.N{i}).(fields.M{j}).uDis]);
                xlabel(ax2,'Dis');
                ylabel(ax2,'Load');
                legend(ax2);
                title(ax2, strcat( fields.M{j}, 'Unknown Input'));
                hold(ax2,'off')
                
                %% Croos_Test (Model Only)
                ax3 = subplot(M,3, 3 + 3*(j-1));
                hold(ax3,'on')
                
                switch(fields.M{j})
                    case {'disR50' 'disR50New'}
                        inputForCross = u50;
                    case {'disR250' 'disR250New'}
                        inputForCross = u250;
                    case {'disR500' 'disR500New'}
                        inputForCross = u500;
                end
                
                for d=1:M
                    if ~strcmp (fields.M{j} , fields.M{d})
                        %Create index vector to sample raw data and match predition data
                        sampleSize = length( inputForCross.dis );
                        tempDis = PmatData.(fields.N{i}).(fields.M{d}).rdis;
                        tempL = length(tempDis);
                        tempIndex = nearest( linspace(1,tempL,sampleSize) );
                        
                        pt = PmatData.(fields.N{i}).(fields.M{d}).rdis (tempIndex ) - ...
                            inputForCross.predictions.(fields.N{i}).(fields.M{d}).model;
                        e = sqrt( mean(pt.^2) );
                        plot(ax3, PmatData.(fields.N{i}).(fields.M{d}).rdis (tempIndex ) , ...
                            PmatData.(fields.N{i}).(fields.M{d}).pload (tempIndex ), ...
                            'DisplayName', strcat( 'ExpData' , fields.M{d} ) );
                        plot(ax3, PmatData.(fields.N{i}).(fields.M{d}).rdis(tempIndex ) , ...
                            inputForCross.predictions.(fields.N{i}).(fields.M{d}).(fields.P{end}),'--', ...
                            'DisplayName', strcat('Model', fields.M{d} ,'@',num2str(e(end),2 )));
                    end
                end
                xlim(ax3,[0 modelFit.(fields.N{i}).(fields.M{j}).uDis]);
                xlabel(ax3,'Dis');
                ylabel(ax3,'Load');
                legend(ax3);
                title(ax3, strcat( fields.M{j}, 'Cross-Input'));
                hold(ax3,'off')
                
            end
            %Enlarge current figure to full size
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            
        end
    otherwise
        disp('Functionality not defined')
end
end

function [NN, NNPerform] = ffnn_single(x,t,T,trainFcn,activationFcn,best_i)
disp("Inside ffnn_single");
%Variables to allocate performance data
NNPerform = zeros(1,T);
NN = cell(1,T);
tr = cell(1,T);
%% This is not required anymore because I am going the data from an
%isolated test.
% Q = size(x, 2);
% Q1 = floor(Q * 0.90);
% Q2 = Q - Q1;
% ind = randperm(Q);
% ind1 = ind(1:Q1);
% ind2 = ind(Q1 + (1:Q2));
% x1 = x(:, ind1);
% t1 = t(:, ind1);
% x2 = x(:, ind2);
% t2 = t(:, ind2);
%%
for j = 1:T
    % Create a Fitting Network
    if nargin < 6
        hiddenLayerSize = j; %Number of neurons
    else
        hiddenLayerSize = [best_i j]; %Number of neurons
    end
    net = feedforwardnet(hiddenLayerSize,trainFcn);
    net.layers{end}.transferFcn = activationFcn;
    % Set up Division of Data for Training, Validation, Testing
    %% When Using TRAINBR Data Division is not Required
    if ~strcmp(trainFcn,'trainbr')
        net.divideParam.trainRatio = 80/100;
        net.divideParam.valRatio = 10/100;
        net.divideParam.testRatio = 10/100;        
    end
    net.userdata = strcat('Neurons:', int2str(j));
    [NN{j},tr{j}] = train(net, x, t);
end
NNPerform = tr;
end

function [neuralnetsbest] = getBestNN(neuralnets)
%Defining number of materials
N = length( fieldnames(neuralnets) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(neuralnets);
neuralnetsbest = neuralnets;

for i=3
    M = length( fieldnames(neuralnets.(fields.N{i})) );
    %Security check when neuralnets is incomplete due to focusing on
    %individual materials
    if M < 2
        continue;
    end
    fields.M = fieldnames(neuralnets.(fields.N{i}));
    for j=6
        neuralnetsbest.(fields.N{i}).(fields.M{j}) = [];
        if strncmp(fields.M(1),'ff',2)
            %New mode: Networks trained with different strain rates(all data compiled)
            [ ~, temp] = ...
                min (neuralnets.(fields.N{i}).(fields.M{j})(1).performance(2:end) );
            index = temp+1;
            %Save network with best performance
            neuralnetsbest.(fields.N{i}).(fields.M{j}).NN = ...
                neuralnets.(fields.N{i}).(fields.M{j})( index).NN;
            neuralnetsbest.(fields.N{i}).(fields.M{j}).index = index;
        else
            %Old mode
            P = length( fieldnames(neuralnets.(fields.N{i}).(fields.M{j}) ) );
            fields.P = fieldnames(neuralnets.(fields.N{i}).(fields.M{j}) );
            for k=1:P
                [ ~, temp] = ...
                    min (neuralnets.(fields.N{i}).(fields.M{j}).(fields.P{k})(1).performance(2:end) );
                index = temp+1;
                
                %Save network with best performance
                neuralnetsbest.(fields.N{i}).(fields.M{j}).(fields.P{k}).NN = ...
                    neuralnets.(fields.N{i}).(fields.M{j}).(fields.P{k})( index).NN;
                neuralnetsbest.(fields.N{i}).(fields.M{j}).(fields.P{k}).index = index;
            end
        end
    end
end
end

function output = nn_eval(input,nn_type_name,NN)
in_s1 = circshift(input,1,2);
in_s1(1) = 0;
in_s2 = circshift(input,2,2);
in_s2(1:2) = 0;
output = 0.*input;
switch nn_type_name
    case {'ff1' 'ff1_m'}
        for i=1:length(input)
            output(i) = NN(input(i));
        end
    case {'ff2' 'ff2_m'}
        for i=1:length(input)
            output(i) = NN( [ input(i); in_s1(i) ] );
        end
    case {'ff3' 'ff3_m'}
        for i=1:length(input)
            output(i) = NN( [ input(i); in_s1(i) ; in_s2(i)] );
        end
    case {'ff4' 'ff4_m'}
        %Initial value
        output(1) = NN( [ input(1); in_s1(1) ;...
            in_s2(1); 0] );
        for i=2:length(input)
            output(i) = NN( [ input(i); in_s1(i) ;...
                in_s2(i); output(i-1)] );
        end
end

end

function output = nn_evalV2(input,disR,nn_type_name,NN)
in_s1 = circshift(input,1,2);
in_s1(1) = 0;
in_s2 = circshift(input,2,2);
in_s2(1:2) = 0;
%Shift second input
in2_s1 = circshift(disR,1,2);
in2_s1(1) = 0;
in2_s2 = circshift(disR,2,2);
in2_s2(1:2) = 0;
output = 0.*input;
switch nn_type_name
    case {'ff1' 'ff1_m'}
        for i=1:length(input)
            output(i) = NN( [ input(i); disR(i) ]);
        end
    case {'ff2' 'ff2_m'}
        for i=1:length(input)
            output(i) = NN( [ input(i); in_s1(i) ; disR(i); in2_s1(i)] );
        end
    case {'ff3' 'ff3_m'}
        for i=1:length(input)
            output(i) = NN( [ input(i); in_s1(i) ; in_s2(i); disR(i); in2_s1(i); in2_s2(i) ] );
        end
    case {'ff4' 'ff4_m'}
        %Initial value
        output(1) = NN( [ input(1); in_s1(1); in_s2(1);...
            disR(1); in2_s1(1); in2_s2(1);...
            0] );
        for i=2:length(input)
            output(i) = NN( [ input(i); in_s1(i); in_s2(i);...
                disR(i); in2_s1(i); in2_s2(i);...
                output(i-1)] );
        end
    case {'ff5' 'ff5_m'}
        %Initial value
        output(1) = NN( [ input(1); disR(1); 0]);
        for i=2:length(input)
            output(i) = NN( [ input(i); disR(i); output(i-1)] );
        end
end

end

function [dis,dt,disR] = defineInput(disR,i_dis,u_dis,size)
%disR input is in mm/min, convert it to mm/sec and return in as output
disR = disR/60; %(mm/min)/60
time = linspace( 0,(u_dis-i_dis)/disR,size);
dt = time(2) - time(1);
dis = time.*disR;
%convert to vector
disR = disR.*ones(1,size);
disR(1) = 0;
end

function [error,P_index,sample_i] = getBestError(config,predictions,PmatData,i,j,fields)

P = length( fieldnames(predictions.(fields.N{i}).(fields.M{j}) ) );
fields.P = fieldnames(predictions.(fields.N{i}).(fields.M{j}) );

%Create index vector to sample raw data and match predition data
sampleSize = length( predictions.(fields.N{i}).(fields.M{j}).(fields.P{1}) );
tempDis = PmatData.(fields.N{i}).(fields.M{j}).rdis;
tempL = length(tempDis);
tempIndex = nearest( linspace(1,tempL,sampleSize) );

for k=1:P
    %Prediction - Target (Error)
    switch(config)
        case 'known'
            pt = predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) - ...
                PmatData.(fields.N{i}).(fields.M{j}).pload (tempIndex);
        case 'unknown'
            pt = predictions.(fields.N{i}).(fields.M{j}).(fields.P{k}) - ...
                predictions.(fields.N{i}).(fields.M{j}).model;
        case 'cross'
            pt = predictions.(fields.N{i}).(fields.M{j}).model - ...
                PmatData.(fields.N{i}).(fields.M{j}).pload (tempIndex);
    end
    %     %Mean Target  - Target (Variation)
    %     mt = mean( PmatData.(fields.N{i}).(fields.M{j}).pload ) - ...
    %         PmatData.(fields.N{i}).(fields.M{j}).pload (tempIndex) ;
    %Root Mean Squared Error
    rmse(k) = sqrt( mean(pt.^2) );
    %             %Relative Squared Error
    %             rse(k) = sum(pt.^2)/sum(mt.^2);
    %             %Mean Absolute Error
    %             mae(k) = mean( abs(pt) );
    %             %Relative Absolute Error
    %             rae(k) = sum( abs(pt) )/sum( abs(mt) );
    %             %Coefficient of Determination (R-Squared)
    %             r2(k) = 1 - ( sum( pt.^2) / sum( mt.^2));
end
%Get Best Error
[~,P_index] = min( rmse(1:(P-1)) );
error = rmse;
sample_i = tempIndex;
end

function getSimulink(net)
%net = nn_best_single.NatR.ff1.NN
gensim(net,-1)
end