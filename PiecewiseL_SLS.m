%Script based on Intron_Tensile.m focusing on the implementation of the
%Piecewise Linearization method on the Standard Linear Solid Model.
%Inputs:
%Smoothed data from Tensile Strength Test stored in TrainingSetProc.mat
%where each material has a dataset for three different deformation ratios:
%50, 250 and 500 mm/min.
%% Initialization
clear all;
clc;
close all;

%force data input (PmatData)
% load('TrainingSetProc.mat'); %Use this path for newly processed data
load('Backup/TrainingSetProc_Sep13.mat'); %Use this path old data
%Defining number of materials
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(PmatData);
%forceing data input (fitParameters)
load('SRelFitParameters.mat');
SRelTests = fieldnames( fitParameters.(fields.N{1})) ;
fields.Q = SRelTests{3};    %Select SRel Test to use L5, L15, L180
%%
%Initializing Core Structure to save Model Fit
modelFit = PmatData;

%% Main loop
for i=1:N
    figure('Name','Fit Accuracy stress vs strain');
    %Secondary for loop accounting for different experimental setups. There
    %are between 3 and 4, material dependent
    M = length( fieldnames( PmatData.(fields.N{i}) ) );
    %Exp setups NAMES: {'disR50';'disR250';'disR500'}
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    %Data variables NAME {'paths';'pforce';'rforce';'rdis';'rtime'}
    P = length( fieldnames( PmatData.(fields.N{i}).(fields.M{1}) ) );
    fields.P = fieldnames( PmatData.(fields.N{i}).(fields.M{1}) );
    for j=1:M
        %%
        %Initializing Core Structure to save Model Fit
        modelFit.( fields.N{i} ).( fields.M{j} ) = [];
        %%
        %Clear variables for each new material
        clear force dis time ke km thao Eo slope strainDepK ;
        %%Definition of ergonomic variables
        force = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{2});
        dis = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{4});
        time = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{5});
        ke = fitParameters.(fields.N{1}).( fields.Q ).SLS.k(1);
        km = fitParameters.(fields.N{1}).( fields.Q ).SLS.k(2);
        thao = fitParameters.(fields.N{1}).( fields.Q ).SLS.thao;
        Eo = fitParameters.(fields.N{1}).( fields.Q ).SLS.Eo;
        [uForce, ult_index] = max(force);
        uDis = dis(ult_index);
        %%
        %%REMOVAL of the equilibrium srping response
        dt = (time(2) - time(1));
        %Extract number of branches, SLSmodel has onlye one branch
        branches = size(thao,2);
        km_force = zeros(branches , length(force));
        for a=2:length(force)
            km_force(branches,a) = (km * (dis(a) - dis(a-1)) + km_force(a-1) ) /...
                (1 + dt/thao);
        end
        real_force = force - km_force;
        %% Decide which part of the curve to analize
        disRange = find(dis < Eo);
        [maxforce , maxforceIndex] = max(force);
%         disRange = 1:maxforceIndex;
        %% Define the variation criteria which will create a new strain
        %segment
        tol = 1; %Percentage
        [divisions, segments] = optimizedStrainSegments(real_force(disRange), dis(disRange), tol);
        
        %% Apply the piecewise linearization
        
        for a=1:length(divisions)-1
            %Linear fit using least mean squares
            fitRange = divisions(a):divisions(a+1);
            %Temporal variables for the stress
            fitReal_force = real_force(fitRange);
            [fit, ~, slope(a) , ~] = linearfit(dis(fitRange),fitReal_force);
            %After the fit the next initial stress point of the line must
            %be equal to the estimated last point of the regression
        end
        %Matrix of coefficients k
        A = zeros(segments,segments);
        b = 1;
        for z=1:segments
            for a=1:b
                A(z,a) = 1;
            end
            b = b + 1;
        end
        PL_ke = A\slope';
        PL_ke = PL_ke';
        %% PIECEWISE LINEARIZATION
        %Initialize variables
        force_PL = zeros(1,length(disRange));
        
        realSumK = PL_ke(1);
        
        d=1;
        for a=3:divisions(end)
            %calculate index of nonlinearK
            if a > divisions(d+1)
                d = d + 1;
                realSumK = realSumK + PL_ke(d);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fit approach  for nonlinear fit and SLS model
            force_PL(a) = ((realSumK + km) * (dis(a) - dis(a-1)) + ...
                realSumK*dis(a)*dt/thao + force_PL(a-1)) / (1 + dt/thao);
        end
        %Subplot of all material
        
        subplot(1,M, j);
        colors = get(gca,'ColorOrder');
        plot(dis(disRange),force(disRange),'Color',colors(1,:));
        hold on
        plot(dis(disRange),force_PL,'--','Color',colors(2,:),'LineWidth',2.5);
        hold off;
        title( strcat( fields.N{i} ,' - ' , fields.M{j}) );
        legend('Data','PL-SLS');
        xlabel('Displacement (mm)');
        ylabel('force (N)');
        grid on;
        ax = gca; % current axes
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.XMinorGrid = 'on';
        ax.YMinorGrid = 'on';
        ax.TickLength = [0.02 0.025];
        %         xlim(parameters(i).xAxisLimits);
        %                 ylim([0 8]);
        %%
        %Good of fitness
        gof(i) = 100*sqrt(mean(power(...
            ( 1/mean(force(disRange)) ) .* (force(disRange) - force_PL),2)));
        
        %%
        %Extraction of key parameters obtained fron PL:
        %segments
        %linearized k's
        %strain segments for linearized k's
        modelFit.( fields.N{i} ).( fields.M{j} ).segments = segments;
        modelFit.( fields.N{i} ).( fields.M{j} ).PL_ke = PL_ke;
        modelFit.( fields.N{i} ).( fields.M{j} ).disIndex = divisions;
        modelFit.( fields.N{i} ).( fields.M{j} ).disSeg = dis(divisions);
        modelFit.( fields.N{i} ).( fields.M{j} ).km = km;
        modelFit.( fields.N{i} ).( fields.M{j} ).thao = thao;
        modelFit.( fields.N{i} ).( fields.M{j} ).Eo = Eo;
        modelFit.( fields.N{i} ).( fields.M{j} ).uDis = uDis;
        modelFit.( fields.N{i} ).( fields.M{j} ).uForce = uForce;
        modelFit.( fields.N{i} ).( fields.M{j} ).ult_index = ult_index;
    end    
end
%%
%Model Fit save into mat file
save('PL.mat','modelFit');

function [segmentsIndex, segments] = optimizedStrainSegments(stress, strain, error)
%PIECEWISE can be implemented by doing consecutives linear regressions and
%analysis the variation in the slope. When a certain tolerance is met, a
%new region/segment


%When the curve slopes varies outside the range +-error. A new strain
%segment is required.error is no in percentage, i.e. it ranges from 0-1
slope = diff(stress)./diff(strain);
%n is the number of founf segments
%j is the index for when a new segment is found, the slope reference changes
%segmentIndex stores the index in which a new segment is found in a form of
%a range hence it starts from 1
n = 1;
j = 1;
segmentsIndex(1) = 1;
for i = 2:length(slope)
    if i < length(slope)
        if slope(i)/slope(j) > (1 + error) || slope(i)/slope(j) < (1 - error)
            n = n + 1;
            segmentsIndex(n) = i; %saves previous and last index to form a range for each segment
            j = i;
        end
    else
        %Last segment
        n = n + 1;
        %diff() reduce the length of stress/strain by 1. In here i+1
        %compensates this
        segmentsIndex(n) = i+1; %saves previous and last index to form a range for each segment
    end
end
segments = n - 1;
% figure
% plot(slope)
end