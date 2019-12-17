%Short version of the ANN_Neuron_Test
%Script to observe the impact of increasing the number of nuerons on a
%feedforward ANN. The matlab tables with the formatted datsets must be
%loaded
% function ANN_Neuron_Test
clear all;
close all;
clc;

load('PiecewiseLinearization.mat'); %For plotting purpouses
load('TrainingSetProc.mat');
known = load('Prediction_KnownDisR_Short.mat');
u50 = load('Prediction_UnknownDisR_50.mat');
u150 = load('Prediction_UnknownDisR_150.mat');
u250 = load('Prediction_UnknownDisR_250.mat');
u400 = load('Prediction_UnknownDisR_400.mat');
u500 = load('Prediction_UnknownDisR_500.mat');


N = length( fieldnames(known.predictions) ); %Value of N is valid for all cases
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(known.predictions);
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