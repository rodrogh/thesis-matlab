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
%Plot visualization
% Change default axes fonts.
set(groot,'DefaultAxesFontName', 'Times New Roman')
set(groot,'DefaultAxesFontSize', 12)
% Change default text fonts.
set(groot,'DefaultTextFontname', 'Times New Roman')
set(groot,'DefaultTextFontSize', 12)
% Extra
set(groot, 'DefaultLegendLocation', 'best')
set(groot, 'defaultLineLineWidth', 2)
%force data input (PmatData)
load('../Experimental work/TrainingSetProcTrainingSetProcStressV4.mat'); %Use this path for newly processed data
% load('Backup/TrainingSetProc_Sep13.mat'); %Use this path old data
%Defining number of materials
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(PmatData);
%forceing data input (fitParameters)
%Version 4 of the data below is pre-Viva
%Version 5 is post-Viva
load('../Experimental work/SRelFitParametersV5.mat');
SRelTests = fieldnames( fitParameters.(fields.N{1})) ;
fields.Q = SRelTests{2};    %Select SRel Test to use L5, L15, L180
%%
%Initializing Core Structure to save Model Fit
modelFit = PmatData;

%% Main loop
%Thesis Plots: Best case is NatRubber with Polyester
%Worst Case is
for i=1:N
    %     figure('Name','Fit Accuracy stress vs strain');
    %Secondary for loop accounting for different experimental setups. There
    %are between 3 and 4, material dependent
    M = length( fieldnames( PmatData.(fields.N{i}) ) );
    %Exp setups NAMES: {'disR50';'disR250';'disR500'}
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    %Data variables NAME {'paths';'pforce';'rforce';'rdis';'rtime'}
    P = length( fieldnames( PmatData.(fields.N{i}).(fields.M{1}) ) );
    fields.P = fieldnames( PmatData.(fields.N{i}).(fields.M{1}) );
    for j=1:M
        %% Only do the strain rate of 500 mm/min for all materials,
        %and 50 mm/min for the SR material
        if strcmp(fields.N{i},"SR" )
            if ~strcmp(fields.M{j},"disR50")
                continue;
            end
        else
            if ~strcmp(fields.M{j},"disR500")
                continue;
            end
        end
        %%
        %Initializing Core Structure to save Model Fit
        modelFit.( fields.N{i} ).( fields.M{j} ) = [];
        %%
        %Clear variables for each new material
        clear slope strainDepK;
        clear force dis time;        
        clear SLS_ke SLS_km SLS_thao Eo;
        clear W_ke W_km W_thao;
        %% Definition of ergonomic variables
        %         force = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{2});
        %         dis = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{4});
        %         time = PmatData.( fields.N{i} ).( fields.M{j} ). (fields.P{5});
        stress = PmatData.( fields.N{i} ).( fields.M{j} ). sstress;
        strain = PmatData.( fields.N{i} ).( fields.M{j} ). pstrain;
        time = PmatData.( fields.N{i} ).( fields.M{j} ). ptime;
        %Stress Relaxation Fit using the Standard Linear Solid Model
        SLS_ke = fitParameters.(fields.N{i}).( fields.Q ).SLS.k(1);
        SLS_km = fitParameters.(fields.N{i}).( fields.Q ).SLS.k(2);
        SLS_thao = fitParameters.(fields.N{i}).( fields.Q ).SLS.thao;
        Eo = fitParameters.(fields.N{i}).( fields.Q ).SLS.Eo;
        %Stress Relaxation Fit using the Generalized Wiechert Model
        W_ke = fitParameters.(fields.N{i}).( fields.Q ).W.k(1);
        W_km = fitParameters.(fields.N{i}).( fields.Q ).W.k(2:end);
        W_thao = fitParameters.(fields.N{i}).( fields.Q ).W.thao;
        %Useful values
        [uStress, uStressIndex] = max(stress);
        uDis = strain(uStressIndex);
        dt = (time(2) - time(1));
        %% Decide which part of the curve to analize: Up to MAX or up to Eo
        %         disRange = find(dis < Eo);
        strainRange = 1:uStressIndex;
        %% Update main variables depending on selected range
        stress = stress(strainRange);
        strain = strain(strainRange);
        %% REMOVAL of the maxwell branches response
        %Obtain cleanded force for SLS model
        SLS_mStress = maxwellStress(strain,SLS_km,SLS_thao,dt);
        SLS_realStress = stress - SLS_mStress;
        %Obtain cleanded force for Generalized Wiechert model
        W_mStress = maxwellStress(strain,W_km,W_thao,dt);
        W_realStress = stress - W_mStress;        
        %% Define SLOPE VARIATION TOLERANCE TO FIT A NEW STRAIN SEGMENT
        tol = 10:10:100; %Percentage
        for t=1:length(tol)
%         for t=1
            [indexOfSegments, valuesOfSegments, NoSegments] = optimizedStrainSegments(SLS_realStress(strainRange), strain(strainRange), tol(t)/100);
            %% PIECEWISE LINEARIZATION ON THE SLS MODEL (OLD)
            %% Apply the piecewise linearization to the obtained segments
            [PL_ke,PL_slope] = piecewiseLin(strain, SLS_realStress, indexOfSegments);            
            %PL_ke contains the values of k1,k2,k3 up to kN
            %PL_slope contains the values of the Strain Dependent Stiffness
            %In other words, k1, k1+k2, k1+k2+k3, k1+k2+k3+...+kN
            stress_PL = zeros(1,length(strainRange));           
            
            for a=2:length(strainRange)
                springIndexes = find( valuesOfSegments < strain(a));
                springDeflection =  strain(a) - valuesOfSegments(springIndexes);
                springResponse = PL_ke(springIndexes).*springDeflection;
                %indexOfSegments has 1 less element than valuesOfSegments
                if isempty(springIndexes)
                    strainDepStiff = PL_slope(1);
                else
                    strainDepStiff = PL_slope( springIndexes(end)); 
                end
                %New approach Jul/01
                stress_PL(a) = (( strainDepStiff + SLS_km) * (strain(a) - strain(a-1)) + ...
                    strainDepStiff*strain(a)*dt/SLS_thao + stress_PL(a-1)) / (1 + dt/SLS_thao);
%                
%                 currentPL = PL_ke_S( dis(:) >= dis(a) );                
%                 force_PL(a) = ((currentPL(1) + SLS_km) * (dis(a) - dis(a-1)) + ...
%                     currentPL(1)*dis(a)*dt/SLS_thao + force_PL(a-1)) / (1 + dt/SLS_thao);
            end
            %% PIECEWISE LINEARIZATION ON THE SLS MODEL (1/Jul)
            %Apply the piecewise linearization to the obtained segments
            [PL_ke,~] = piecewiseLin(strain, SLS_realStress, indexOfSegments);           
            [SDSStress] = SDStiffnessStress(strain,PL_ke,valuesOfSegments);
            stress_SLS_PL = SDSStress + SLS_mStress;
            %% PIECEWISE LINEARIZATION ON THE WIECHERT MODEL (1/Jul)
            %Apply the piecewise linearization to the obtained segments
            [PL_ke,~] = piecewiseLin(strain, W_realStress, indexOfSegments);           
            [SDSStress] = SDStiffnessStress(strain,PL_ke,valuesOfSegments);
            stress_W_PL = SDSStress + W_mStress;
            %Subplot of all material
            colors = get(gca,'ColorOrder');
            switch( fields.M{j} )
                case {'disR50','disR50New'}
                    test_title = "50 mm/min";
                case {'disR250','disR250New'}
                    test_title = "250 mm/min";
                case {'disR500','disR500New'}
                    test_title = "500 mm/min";
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %This plot are the standard
            %         subplot(1,M, j);
%             figure
% %             colors = get(gca,'ColorOrder');
%             plot(strain,stress*1e-6,'Color',colors(6,:),'LineWidth',1.5);
%             hold on
%             plot(strain,force_SLS_PL*1e-6,':','Color',colors(7,:),'LineWidth',1.5);
%             plot(strain,force_W_PL*1e-6,':','Color',colors(2,:),'LineWidth',1.5);
%             hold off;
%             title( strcat( fields.N{i} ," - " , test_title) );
%             legend('Experimental Data','PL-SLS Fit','PL-W Fit');
%             xlabel('Strain');
%             ylabel('Stress (MPa)');
%             grid on;
%             ax = gca; % current axes
%             ax.Box = 'off';
%             ax.XMinorTick = 'on';
%             ax.YMinorTick = 'on';
%             ax.XMinorGrid = 'on';
%             ax.YMinorGrid = 'on';
%             ax.TickLength = [0.02 0.025];
% %                     xlim(parameters(i).xAxisLimits);
% %                             ylim([0 8]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Good of fitness (Normalized RMSE)
            gof = 100*sqrt( mean( (stress - stress_PL).^2 )./ mean( stress.^2 ) );
            gofNMAD = 100*( mean( abs( stress(strainRange) - stress_PL) ) )/ mean(stress(strainRange));
            %New approaches (Jul/1)
            gof_SLS = 100*sqrt( mean( (stress - stress_SLS_PL).^2 )./ mean( stress.^2 ) );
%             gof_SLS = sqrt( mean( (stress - stress_SLS_PL).^2 ) );
            gof_W = 100*sqrt( mean( (stress - stress_W_PL).^2 )./ mean( stress.^2 ) );            
%             gof_W = sqrt( mean( (stress - stress_W_PL).^2 ) );
            modelFit.( fields.N{i} ).( fields.M{j} )(t).gofSLS = gof_SLS;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).gofW = gof_W;
            %% Extraction of key parameters obtained fron PL:
            %segments
            %linearized k's
            %strain segments for linearized k's
            modelFit.( fields.N{i} ).( fields.M{j} )(t).tol = tol(t);
            modelFit.( fields.N{i} ).( fields.M{j} )(t).gof = gof;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).gofNMAD = gofNMAD;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).segments = NoSegments;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).PL_ke = PL_ke;
%             modelFit.( fields.N{i} ).( fields.M{j} )(t).PL_ke_S = PL_ke_S;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).PL_ke_dis = strain; %This variable relates the PL_ke_S to the fitted distance
            modelFit.( fields.N{i} ).( fields.M{j} )(t).disIndex = indexOfSegments;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).disSeg = strain(indexOfSegments);
            modelFit.( fields.N{i} ).( fields.M{j} )(t).km = SLS_km;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).thao = SLS_thao;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).Eo = Eo;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).uDis = uDis;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).uForce = uStress;
            modelFit.( fields.N{i} ).( fields.M{j} )(t).uForceIndex = uStressIndex;
        end
    end
end
%% Plot tradeoff
plotPLTradeOff(modelFit)
%% Save Model Fit
%save('PL.mat','modelFit');

function plotPLTradeOff(modelFit)
%% Plots: Trade-off between complexity and accuracy
%Defining number of materials
N = length( fieldnames(modelFit) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(modelFit);
for i=1:N
    %Secondary for loop accounting for different experimental setups. There
    %are between 3 and 4, material dependent
    M = length( fieldnames( modelFit.(fields.N{i}) ) );
    %Exp setups NAMES: {'disR50';'disR250';'disR500'}
    fields.M = fieldnames( modelFit.(fields.N{i}) );    
    for j=1:M
        %% Only do the strain rate of 500 mm/min for all materials,
        %and 50 mm/min for the SR material
        if strcmp(fields.N{i},"SR" )
            if ~strcmp(fields.M{j},"disR50")
                continue;
            end
        else
            if ~strcmp(fields.M{j},"disR500")
                continue;
            end
        end
        P = length(modelFit.( fields.N{i} ).( fields.M{j} ));
        gofOld = zeros(1,P);
        gofSLS = zeros(1,P);
        gofW = zeros(1,P);
        NoSegments = zeros(1,P);
        tol = zeros(1,P);
        for k=1:P
            gofOld(k) = modelFit.( fields.N{i} ).( fields.M{j} )(k).gof;
            gofSLS(k) = modelFit.( fields.N{i} ).( fields.M{j} )(k).gofSLS;
            gofW(k) = modelFit.( fields.N{i} ).( fields.M{j} )(k).gofW;
            NoSegments(k) = modelFit.( fields.N{i} ).( fields.M{j} )(k).segments;
            tol(k) = modelFit.( fields.N{i} ).( fields.M{j} )(k).tol;
        end
    end
    figure
    yyaxis left
    bar(tol,NoSegments)
    
    yyaxis right
    semilogy(tol,gofSLS,tol,gofW);
    
    legend;
end

end

function [mStress] = maxwellStress(strain,k,thao,dt)
%This finction returns a vector with all the values of the stress response
%in the maxwell branches depending on the applied strain.
%The k and thao parameters are from the maxwell branches, as many as
%available
N = length(strain);
M = length(thao);
stress = zeros(M,N);
for i=2:N %Stress when strain=0 is 0
    for j =1:M
        %Matrix of individual stress response per spring
        stress(j,i) = (k(j)*(strain(i) - strain(i-1)) + stress(j,i-1) ) /...
                   (1 + dt/thao(j));
    end
end
mStress = sum(stress,1);
end

function [PL_ke,PL_slope] = piecewiseLin(strain, real_stress, indexOfSegments)
%% Apply the piecewise linearization            
NoSegments = length(indexOfSegments)-1;
PL_slope = zeros(1,NoSegments);
for a=1:NoSegments
    %Linear fit using least mean squares
    fitRange = indexOfSegments(a):indexOfSegments(a+1);
    [~, ~, PL_slope(a) , ~] = linearfit(strain(fitRange),real_stress(fitRange));
    %After the fit the next initial stress point of the line must
    %be equal to the estimated last point of the regression
end
%Matrix of coefficients k
A = zeros(NoSegments,NoSegments);
b = 1;
for z=1:NoSegments
    for a=1:b
        A(z,a) = 1;
    end
    b = b + 1;
end
PL_ke = (A\PL_slope')';
end

function [totalStress] = SDStiffnessStress(strain,PL_ke,valuesOfSegments)
%Stress reponse of the Strain Dependent Stiffness
N = length(strain);
totalStress = strain.*0;
for i = 1:N
    indexes = find( valuesOfSegments < strain(i));
    if isempty(indexes)
        indexes = 1;
    end
    indivStrain =  strain(i) - valuesOfSegments(indexes);
    indivStress = PL_ke(indexes).*indivStrain;
    totalStress(i) = sum(indivStress,2);
end
end

function [indexOfSegments,valuesOfSegments,noSegments] = optimizedStrainSegments(stress,strain,error)
%Inputs:
%stress - 1xN vector
%strain - 1xN vector
%error - scalar values from (0.01:1)
%Outputs:

%Auxiliary variables
N = length(stress);
toleranceMet = zeros(1,N); 
%Nonzero values indicate indexes of obtained strain segments
toleranceMet(1) = 1;
%Calculate stress-strain curve slope
slope = diff(stress)./diff(strain); %vector 1x(N-1)
%Consider Adding Smooth to the slope to avoid cluttering the identified
%stratin segments at the beginning of the curve

%Define reference value for searching algorithm
slopeRef = slope(1); %scalar

%Iterative search 
for i=1:N-1 %vector 1x(N-1)
    %Calculate slope changes with respect to the reference
    slopeChange = abs( (slopeRef-slope(i) ) / slopeRef);
    %Find changes greater than error
    if slopeChange >= error
        slopeRef = slope(i); %update reference
        toleranceMet(i) = 1;
    end
end
%Nonzero values indicate indexes of obtained strain segments
indexOfSegments = find(toleranceMet);
%In the event that the iterative search found no viable segments
if length(indexOfSegments) == 1
    segments = [strain(1),strain(end)];
else
    segments = [strain(indexOfSegments),strain(end)];
end

indexOfSegments = [indexOfSegments,N]; %Add value of last index lost during diff
valuesOfSegments = segments;
noSegments = length(segments)-1;
end

function [segmentsIndex, segments] = OldoptimizedStrainSegments(stress, strain, error)
%PIECEWISE can be implemented by doing consecutives linear regressions and
%analysis the variation in the slope. When a certain tolerance is met, a
%new region/segment


%When the curve slopes varies outside the range +-error. A new strain
%segment is required.error is no in percentage, i.e. it ranges from 0-1
slope = diff(stress)./diff(strain);
% %n is the number of founf segments
% %j is the index for when a new segment is found, the slope reference changes
% %segmentIndex stores the index in which a new segment is found in a form of
% %a range hence it starts from 1
% n = 1;
% j = 1;
% segmentsIndex(1) = 1;
% for i = 2:length(slope)
%     if i < length(slope)
%         if slope(i)/slope(j) > (1 + error) || slope(i)/slope(j) < (1 - error)
%             n = n + 1;
%             segmentsIndex(n) = i; %saves previous and last index to form a range for each segment
%             j = i;
%         end
%     else
%         %Last segment
%         n = n + 1;
%         %diff() reduce the length of stress/strain by 1. In here i+1
%         %compensates this
%         segmentsIndex(n) = i+1; %saves previous and last index to form a range for each segment
%     end
% end
% segments = n - 1;
%% New Approach
%Numerical differentiation of the slope substract 1 to the dimension of the
%used variables. This needs to be compensated
refSlope = slope(1);
segmentsIndex = 1;
while true    
    newSlope = refSlope./slope;
    slopePosChange = find( newSlope > (1 + error));
    slopeNegChange = find( newSlope < (1 - error));
    %Find which change, either positive or negative, happens first by
    %calculating the min of the first index of both variables
    minPosChange = min(slopePosChange(segmentsIndex(end):end ),[],'omitnan' );
    minNegChange = min(slopeNegChange(segmentsIndex(end):end ),[],'omitnan' );
    if ~isempty(minPosChange) || ~isempty(minNegChange)        
        newSegmentIndex = min( [minPosChange,minNegChange],[],'omitnan');
        segmentsIndex = [segmentsIndex,newSegmentIndex];
        refSlope = slope(newSegmentIndex);
    else
        if length(segmentsIndex) == 1
            segmentsIndex = [1,length(slope)+1];
        else
            if segmentsIndex(end) < (length(slope)+1)
                segmentsIndex = [segmentsIndex,length(slope)+1];
            end
        end        
        segments = length(segmentsIndex)-1;
        break;
    end    
end

end