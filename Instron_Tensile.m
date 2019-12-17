%Postprocessing of Instron Tensile Experiment data
%The first important feature to obtain is the Yield Strength of the
%material in order to stablish a linear region where the material should be
%when performing time-dependant experiments
clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The data is distributed in six columns as follows:
%1 - Time
%2 - Extension
%3 - Load
%4 - Tensile Strain
%5 - Tensile Stress
%6 - Corrected position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Font parameters for plotting
lineSize = 3;
fontSize = 18;
font = 'Gill Sans MT';
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lineSize);   % set the default line width to lw
set(0,'defaultAxesFontName',font);   % set the default line width to lw
set(0,'defaultAxesFontSize',fontSize);   % set the default line width to lw
set(0,'defaultTextFontName',font);   % set the default line width to lw
set(0,'defaultTextFontSize',fontSize);   % set the default line width to lw
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter for data extraction. Acces each element with filePaths(n).paths
N = 6;  %Number of experiments
lo = 33; %initial length in mm since experimental data is in mm as well
parameters = struct('name',...
    {'Ethylene Polypropylene Rubber';
    'Fluorocarbon Rubber';
    'Natural Rubber';
    'Nitrile Rubber';
    'Polyethylene Rubber';
    'Silicone Rubber';},...
    'paths',...
    {'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_'},...
    'shortName',...
    {'EPR';
    'FR';
    'NatR';
    'NR';
    'PR';
    'SR'},...
    'NNFunction',...
    {'NeuralNetwork_EPR';
    'NeuralNetwork_FR';
    'NeuralNetwork_NatR';
    'NeuralNetwork_NR';
    'NeuralNetwork_PR6';
    'NeuralNetwork_SR'},...
    'NN_P_Function',...
    {'NeuralNetwork_P_EPR';
    'NeuralNetwork_P_FR';
    'NeuralNetwork_P_NatR';
    'NeuralNetwork_P_NR';
    'NeuralNetwork_P_PR6';
    'NeuralNetwork_P_SR'},...
    'NN_AppP_Function',...
    {'NeuralNetwork_AppP_EPR';
    'NeuralNetwork_AppP_FR';
    'NeuralNetwork_AppP_NatR';
    'NeuralNetwork_AppP_NR';
    'NeuralNetwork_AppP_PR';
    'NeuralNetwork_AppP_SR'},...
    'toeRegion',...
    {1:15;  1:5;    1:28;   1:6;    1:5;    1:8},...
    'elRegion',...
    {5:61;  5:61;   28:841; 6:73;   1:60;  8:721},...
    'toe',...
    {9;  8;   3; 8;   4;  9},...
    'el',...
    {20;  20;   20; 20;   15;  20},...
    'area',... %in meters since load is given in Newtons
    {(6e-3*1.5e-3);    (6e-3*1.5e-3);    (6e-3*1.5e-3);    (6e-3*1.5e-3);  (6e-3*6e-3);  (6e-3*1.5e-3)},...
    'strainRate',... %in meters per second
    {500*(1e-3/60); 500*(1e-3/60); 50*(1e-3/60); 500*(1e-3/60); 500*(1e-3/60); 50*(1e-3/60)},...
    'strainRelax',...
    {5e-3/lo; 5e-3/lo; 7e-3/lo; 6e-3/lo; 3e-3/lo; 6e-3/lo},...
    'xAxisLimits',...
    {[0 1000]; [0 500]; [0 150]; [0 550]; [0 240]; [0 750]},...
    'yAxisLimits',...
    {[0 80]; [0 45]; [0 35]; [0 40]; [0 12]; [0 65]});
%{20e-3/lo; 10e-3/lo; 6e-3/lo; 5e-3/lo; 4e-3/lo; 15e-3/lo});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Elastic Limit of all materials
%EPR - Elastic region before 5mm elongation
%FR - Elastic region before 5mm elongation
%NatR - Elastic region before 7mm elongation
%NR - Elastic region before 6mm elongation
%PR3 - Elastic region before 3mm elongation
%PR6 - Elastic region before 3mm elongation
%SR - Elastic region before 6mm elongation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Figures created
fignStress = figure('Name','Normalized stress vs Normalized strain');
% figStiff = figure('Name','Stiffness vs Strain')
figComparison = figure('Name','Fit Accuracy stress vs strain');
figTrousers = figure('Name','Tendon vs Rubbers');
figIMMX = figure('Name','Fit Accuracy stress vs strain');
figgof = figure('Name','Fit Accuracy stress vs strain');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Row and Column where data is being held
iRow = 2;
iCol = 0;
files = 5; % Number of files to read
expData = struct('data',{},'mean',{},'strain',{},'stress',{},'YMtoe',{},'YM',{},...
    'uStrain',{},'uStress',{},'uDis',{},'uLoad',{},'load',{},'dis',{},...
    'stiffToe',{},'stiffEl',{});

MatPar = struct('k',{},'thao',{},'PL_k',{},'strain_PL_k',{});

%Load data processed from readInstronTable.m
%To create new dataset and then load it. Comment this line an uncommentthe
%last line of the code. Then reverse.
% load('readInstronData.mat');
%Load data from Stress-Rel Fit using Maxwell model equations. The file
%contains k and eta. Low or High refers to the maxmum number of branches
%used in the fitting
% load('FitDataSRLow.mat');
for i=1:N
    [expData(i).load, ~ , expData(i).dis, expData(i).time] = readInstronTable(files,iRow,iCol,parameters(i).paths,'TensileStrength');
    continue;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Obtaining Yield Strength. The linear elastic region of each material is
    %selected by visual observation. The slope of this portion of the data is
    % obtained using linear regression. The 0.2% of elongation method is applied
    % to obtain the Yield Stress.
    %Create more ergonomic variables for each parameter
    time = expData(i).time;
    dis = expData(i).dis;
    loadn = expData(i).load;
    expData(i).strain = expData(i).dis./lo;
    expData(i).stress = expData(i).load./(parameters(i).area);
    %Ergonomic varuables
    strain = expData(i).strain;
    stress = expData(i).stress;
    %Ultimate mechanical properties
    [expData(i).uStress, maxStressIndex] = max(stress);
    expData(i).uStrain = strain(maxStressIndex);
    [expData(i).uLoad, maxLoadIndex] = max(expData(i).load);
    expData(i).uDis = expData(i).dis(maxLoadIndex);    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The nonlinear stiffness of the rubbers will be linearized with a piecewise
    %algorithm. Which means that N diferent stiffness will be extracted from
    %each material at N different strains. The assumption is that the stiffness
    %ki is a contribution of previous stiffnesses.
    %Piecewise fit of the nonlinear stiffness of the system. Knowing all
    %the coefficients k, the viscous element coefficient can be known
    %Definition of strain step. The Entirity of the tensile strength data
    %must be divided in N parts to be fitted.
    %Fit data up to the strain used in StressRel test
    %     fitNumbers = 51;
    %for fitNumbers = 3:3
    elements = 0;
    for e=1:elements+1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Decide which part of the curve to analize
%         strainRange = find(strain < Eo(i));
        strainRange = 11:maxStressIndex;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Accounting for the stress in the branches have proven benefitial for
        %the SLS model. Using Wiechert model fit in Srel with n branches
        dt = (time(2) - time(1));
        branches = size(eta,2);
        stressBSum = zeros(size(stress));
        stressBTemp = zeros(length(stress),branches);
        for m=2:length(stress)
            for n=1:Branches(i)
                stressBTemp(m,n) = (k(i,n+1) * (strain(m) - strain(m-1)) + stressBTemp(m-1,n) ) /...
                    (1 + dt/thao(i,n));
                stressBSum(m) = stressBSum(m) + stressBTemp(m,n);
            end
        end
        realStress = stress - stressBSum;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Accounting for the stress in the branches have proven benefitial for
        %the SLS model. Using SLS fit in Srel with 1 branch
        dt = (time(2) - time(1));
        branches = size(etaSLS,2);
        stressBSum = zeros(size(stress));
        stressBTemp = zeros(length(stress),branches);
        for m=2:length(stress)
            stressBTemp(m) = (kSLS(i,2) * (strain(m) - strain(m-1)) + stressBTemp(m-1) ) /...
                (1 + dt/thaoSLS(i));
            stressBSum(m) = stressBSum(m) + stressBTemp(m);
        end
        SLSStress = stress - stressBSum;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The new Wiechert model with 2 term will only account for the stress
        %offset of those two branches, as follows
        dt = time(2) - time(1);
        branches = size(etaW,2);
        stressBSum = zeros(size(stress));
        stressBTemp = zeros(length(stress),branches);
        for m=2:length(stress)
            for n=1:branches
                stressBTemp(m,n) = (kW(i,n+1) * (strain(m) - strain(m-1)) + stressBTemp(m-1,n) ) /...
                    (1 + dt/thaoW(i,n));
                stressBSum(m) = stressBSum(m) + stressBTemp(m,n);
            end
        end
        WStress = stress - stressBSum;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Define the variation criteria which will create a new strain
        %segment
        initTol = 10;
        endTol = 0;
        error(e) = (initTol + (e-1)*(endTol-initTol)/elements)/100; %Error tolerance is %
        [divisions, segments(i,e)] = optimizedStrainSegments(realStress(strainRange), strain(strainRange), 0.1);
        %         divisions = round(linspace(1,strainRange(end),fitNumbers));
        %         figure
        %         plot(strain*100,1e-6*stress);
        %         hold on
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Apply the piecewise linearization
        clear slope;
        clear slope1;
        clear slopeW;
        clear slopeSLS;
        for j=1:length(divisions)-1
            %Linear fit using least mean squares
            fitRange = divisions(j):divisions(j+1);
            %Temporal variables for the stress
            fitStress = stress(fitRange);
            fitRealStress = realStress(fitRange);
            fitWStress = WStress(fitRange);
            fitSLSStress = SLSStress(fitRange);
            [fit, ~, slope(j),~] = linearfit(strain(fitRange),fitStress);
            [fitReal, ~, slope1(j),~] = linearfit(strain(fitRange),fitRealStress);
            [fitW, ~, slopeW(j),~] = linearfit(strain(fitRange),fitWStress);
            [fitSLS, ~, slopeSLS(j),~] = linearfit(strain(fitRange),fitSLSStress);
            %             plot(100*strain(fitRange),1e-6*fit);
            %After the fit the next initial stress point of the line must
            %be equal to the estimated last point of the regression
            
        end
        %Matrix of coefficients k
        %         A = zeros(fitNumbers-1,fitNumbers-1);
        A = zeros(segments(i,e),segments(i,e));
        b = 1;
        %         for j=1:fitNumbers-1
        for j=1:segments(i,e)
            for a=1:b
                A(j,a) = 1;
            end
            b = b + 1;
        end
        nonlinearK = A\slope';
        realNonlinearK = A\slope1';
        WNonlinearK = A\slopeW';
        SLSNonlinearK = A\slopeSLS';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Use the fitting approach in your notes to obtain the stress at a given
        %moment and compare it with the experimental data. This is a time
        %dependant function, hence the time vector is required
        %Initial conditions for stress at the first two instants
        d = 1;
        dt = time(2) - time(1);
        %                 dt = 1;
        branches = size(eta,2);
        %Different fitting approaches
        %Fit 1 = combined nonlinear + wiechert model fit with the equation
        %solution for constant strain
        %Fit 2 = combined nonlinear + wiechert model fit with general equation
        %Fit 3 = wiechert model fit
        %Fit 4 = nonlinear linearization with Paper equation which relates the
        %sitffness K in SLS model with the strain rate
        stressFit_1 = zeros(length(strainRange),1);
        stressFit_2 = zeros(length(strainRange),1);
        stressFit_4 = zeros(length(strainRange),1);
        stressFit_4W = zeros(length(strainRange),1);
        stressFit_4SLS = zeros(length(strainRange),1);
        stressFit_5 = zeros(length(strainRange),1);
        stressB_1 = zeros(length(strainRange),1);
        stressB_2 = zeros(length(strainRange),branches);
        stressNN = zeros(length(strainRange),1);
        stressNN_multi = zeros(length(strainRange),1);
        stressRNN_RI = zeros(length(strainRange),1);
        stressRNN_RD = zeros(length(strainRange),1);
        %First k of the k vector obtained from fit is ke
        ke = k(i,1);
        
        sumK = nonlinearK(1);
        realSumK = realNonlinearK(1);
        WSumK = WNonlinearK(1);
        SLSSumK = SLSNonlinearK(1);
        for j=3:divisions(end)
            %calculate index of nonlinearK
            if j > divisions(d+1)
                d = d + 1;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Correction of the heavise function to improve response
                %Approximate the 1/2 heaviside equivalent of the stress
                %response when stiffness change
                %Wiechert - Equation for unkown strain input
                stressFit_2(j-1) = strain(j-1)*(realSumK + 0.5*realNonlinearK(d));
                for n=1:Branches(i)
                    stressB_2(j-1,n) = (k(i,n+1) * (strain(j-1) - strain(j-2)) + stressB_2(j-2,n) ) /...
                        (1 + dt/thao(i,n));
                    stressFit_2(j-1) = stressFit_2(j-1) + stressB_2(j-1,n);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                sumK = sumK + nonlinearK(d);
                realSumK = realSumK + realNonlinearK(d);
                WSumK = WSumK + WNonlinearK(d);
                SLSSumK = SLSSumK + SLSNonlinearK(d);
            end
            %Wiechert - Equation for unkown strain input
            stressFit_2(j) = strain(j)*realSumK;
            for n=1:Branches(i)
                stressB_2(j,n) = (k(i,n+1) * (strain(j) - strain(j-1)) + stressB_2(j-1,n) ) /...
                    (1 + dt/thao(i,n));
                stressFit_2(j) = stressFit_2(j) + stressB_2(j,n);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fit approach  for nonlinear fit and SLS model using Wiechert
            %branches
            stressFit_4(j) = ((realSumK + kSLS(i,2)) * (strain(j) - strain(j-1)) + ...
                realSumK*strain(j)*dt/thaoSLS(i) + stressFit_4(j-1)) / (1 + dt/thaoSLS(i));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fit approach  for nonlinear fit and SLS model using SLS
            %branches
            stressFit_4SLS(j) = ((SLSSumK + kSLS(i,2)) * (strain(j) - strain(j-1)) + ...
                SLSSumK*strain(j)*dt/thaoSLS(i) + stressFit_4SLS(j-1)) / (1 + dt/thaoSLS(i));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fit approach  for nonlinear fit and SLS model using 2
            %branches
            stressFit_4W(j) = ((WSumK + kSLS(i,2)) * (strain(j) - strain(j-1)) + ...
                WSumK*strain(j)*dt/thaoSLS(i) + stressFit_4SLS(j-1)) / (1 + dt/thaoSLS(i));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Fit approach  for nonlinear fit and wiecher model with 2 branches
            %             strainDepK = sumK;
            strainDepK = realSumK;
            %                         strainDepK = WSumK;
            %                         strainDepK = kW(i,1);
            
            %Defining constants for Wiechert model 2 terms
            %         a2 = thao(i,1)*thao(i,branches);
            %         a1 = thao(i,1) + thao(i,branches);
            %         a0 = 1;
            %         alpha = a2/(dt^2);
            %         b2 = thao(i,1)*thao(i,branches) * ( sumK + k(i,2) + k(i,branches));
            %         b1 = sumK * (thao(i,1) + thao(i,branches)) + k(i,2)*thao(i,1) + k(i,branches)*thao(i,branches);
            %         b0 = sumK;
            %         beta = b2/(dt^2);
            %             a2 = thaoW(i,1)*thaoW(i,2);
            %             a1 = thaoW(i,1) + thaoW(i,2);
            %             a0 = 1;
            %             alpha = a2/(dt^2);
            %             b2 = thaoW(i,1)*thaoW(i,2) * ( strainDepK + kW(i,2) + kW(i,3));
            %             b1 = strainDepK * (thaoW(i,1) + thaoW(i,2)) + kW(i,2)*thaoW(i,1) + kW(i,3)*thaoW(i,2);
            %             b0 = strainDepK;
            %             beta = b2/(dt^2);
            %New constants
            %             a2 = thaoW(i,1)*thaoW(i,2);
            %             a1 = kW(i,2)*thaoW(i,2) + kW(i,3)*thaoW(i,1);
            %             a0 = kW(i,2)*kW(i,3);
            %             b2 = thaoW(i,1)*thaoW(i,2) * ( strainDepK + kW(i,2) + kW(i,3));
            %             b1 = kW(i,2)*kW(i,3)*(thaoW(i,1) + thaoW(i,2)) + strainDepK * (kW(i,2)*thaoW(i,2) + kW(i,3)*thaoW(i,1));
            %             b0 = strainDepK*kW(i,2)*kW(i,3);
            %             c2 = a2/(dt*dt);
            %             c1 = 2*a2/(dt*dt) + a1/dt;
            %             c0 = a2/(dt*dt) + a1/dt + a0;
            %             %Initializing finite differences equations Wiechert 2 branches
            %             %             stressFit_5(2) = (stressFit_5(j-1)*c1...
            %             %                     + ((strain(j) - 2*strain(j-1))/(dt*dt)) * b2 ...
            %             %                     + ((strain(j) - strain(j-1))/dt) * b1 ...
            %             %                     + strain(j)*b0)...
            %             %                     / c0;
            %
            %             %             if j > 2
            %             straind2 = (strain(j) - 2*strain(j-1) + strain(j-2))/(dt*dt);
            %             straind1 = (strain(j) - strain(j-1))/dt;
            %
            %             stressFit_5(j) = (stressFit_5(j-1)*c1...
            %                 - stressFit_5(j-2)*c2...
            %                 + straind2*b2 + straind1*b1 + strain(j)*b0)...
            %                 / c0;
            %             %             end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Prepare data to test Neural Networks
        strainRate = diff(strain(strainRange))./diff(time(strainRange));
        strainRate = [strainRate ; strainRate(end)];
        %The network receives distances and provide loads as output
        neuralInput = lo.*[strain(strainRange) strainRate]';
        %Prepare data for shifted networks
        temp = neuralInput';
        s1 = circshift(temp,1,1); %Shift both dimensions by 1
        s1(1,:) = 0;
        s2 = circshift(temp,2,1); %Shift both dimensions by 2
        s2(1:2,:) = 0;
        %Save shifted values in order, first displacement and its t,t-1,t-2
        inputShift = [temp(:,1),s1(:,1),s2(:,1),...
            temp(:,2),s1(:,2),s2(:,2)]';
        NN_fh = str2func( strcat('NN1H_P_',parameters(i).shortName) );
        NN_multi_fh = str2func( strcat('NN2H_P_',parameters(i).shortName) );
        RNN_RI_fh = str2func( strcat('RNN_RI_1S_O_P_',parameters(i).shortName) );
        RNN_RD_fh = str2func( strcat('RNN_RD_1S_O_P_',parameters(i).shortName) );
        
        stressNN = NN_fh(neuralInput)/parameters(i).area;
        stressNN_multi = NN_multi_fh(neuralInput)/parameters(i).area;
        initO_RI = 0;
        initO_RD = 0;
        for x=1:size(neuralInput(strainRange),2)
            formatI = [ inputShift( 1 , x) ;...
                inputShift( 2 , x) ;...
                initO_RI];
            stressRNN_RI(x) = RNN_RI_fh( formatI )/parameters(i).area;
            initO_RI = stressRNN_RI(x);
            clear formatI
            formatI = [ inputShift( 1 , x) ;...
                inputShift( 2 , x) ;...
                inputShift( 4 , x) ;...
                inputShift( 5 , x) ;...
                initO_RD];
            stressRNN_RD(x) = RNN_RD_fh( formatI )/parameters(i).area;
            initO_RD = stressRNN_RD(x);
            clear formatI
        end
        %Subplot of all material
        figure(figComparison)
        subplot(2,3,i);
        %         figure
        colors = get(gca,'ColorOrder');
        plot(strain(strainRange)*100,1e-6*stress(strainRange),'Color',colors(6,:));
        hold on
        %         plot(strain(strainRange)*100,1e-6*stressFit_2,'Color',colors(5,:));
        plot(strain(strainRange)*100,1e-6*stressFit_4,'--','Color',colors(7,:),'LineWidth',2.5);
        plot(strain(strainRange)*100,1e-6*stressFit_4SLS)
        plot(strain(strainRange)*100,1e-6*stressFit_4W,'--')
%         plot(strain(strainRange)*100,1e-6*stressNN);
%         plot(strain(strainRange)*100,1e-6*stressNN_multi);
%         plot(strain(strainRange)*100,1e-6*stressRNN_RI);
%         plot(strain(strainRange)*100,1e-6*stressRNN_RD);
        %                 plot(strain(strainRange)*100,1e-6*stressFit_5);
        hold off
        sTitle = strcat(parameters(i).name,'');
        legend('Data','PL-SLS','PL-SLS-1','PL-SLS-2');
%         legend('Data','PL-SLS','NN','NN_multi','RNN_RI','RNN_RD');
        title(sTitle);
        xlabel('Strain (%)');
        ylabel('Stress (MPa)');
        grid on;
        ax = gca; % current axes
        ax.XMinorTick = 'on';
        ax.YMinorTick = 'on';
        ax.XMinorGrid = 'on';
        ax.YMinorGrid = 'on';
        ax.TickLength = [0.02 0.025];
        %         xlim(parameters(i).xAxisLimits);
        %                 ylim([0 8]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Good of fitness
        gofW(i,e) = 100*sqrt(mean(power(...
            (1/mean(stress(strainRange)))*(stress(strainRange) - stressFit_2),2))); %Wiechert
        gofSLS(i,e) = 100*sqrt(mean(power(...
            (1/mean(stress(strainRange)))*(stress(strainRange) - stressFit_4),2))); %SLS
        gofSLS_1(i,e) = 100*sqrt(mean(power(...
            (1/mean(stress(strainRange)))*(stress(strainRange) - stressFit_4SLS),2))); %SLS
        gofSLS_2(i,e) = 100*sqrt(mean(power(...
            (1/mean(stress(strainRange)))*(stress(strainRange) - stressFit_4W),2))); %SLS
%         gofNN(i,e) = 100*sqrt(mean(power(...
%             (1/mean(stress(strainRange)))*(stress(strainRange) - stressNN'),2)));
%         gofNN_multi(i,e) = 100*sqrt(mean(power(...
%             (1/mean(stress(strainRange)))*(stress(strainRange) - stressNN_multi'),2)));
%         gofRNN_RI(i,e) = 100*sqrt(mean(power(...
%             (1/mean(stress(strainRange)))*(stress(strainRange) - stressRNN_RI),2)));
%         gofRNN_RD(i,e) = 100*sqrt(mean(power(...
%             (1/mean(stress(strainRange)))*(stress(strainRange) - stressRNN_RD),2)));
    end
    MatPar(i).k = kSLS(i,2);
    MatPar(i).thao = thao(i);
    MatPar(i).strain_PL_k = strain(divisions);
    MatPar(i).PL_k = slopeSLS';
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %         Plot of the GOF per material along different number of segments
    %This plot only works when different number of segments are studied
    %     figure(figgof)
    %     subplot(2,3,i);
    % %     figure
    %     colors = get(gca,'ColorOrder');
    %     yyaxis left
    %     bar(error*100,segments(i,:),'FaceColor',colors(6,:))
    %     ylabel('Strain Segments');
    %
    %     yyaxis right
    %     plot(error*100,gofSLS(i,:),'-','Color',colors(7,:));
    %     hold on
    %     plot(error*100,gofW(i,:),'-','Color',colors(5,:));
    %     %                 plot(gofW2(i,:));
    %     hold off
    %     legend('Strain Segments','PL - SLS','PL - Wiechert');
    %     sTitle = parameters(i).name;
    %     title(sTitle);
    %     ylabel('Relative RMSE (%)');
    %     xlabel('Slope Variation Tolerance (%)');
    %     set(gca, 'XDir','reverse')
    %     xlim([0 105]);
    %     grid on;
    %     ax = gca; % current axes
    %     ax.XMinorTick = 'on';
    %     ax.YMinorTick = 'on';
    %     ax.XMinorGrid = 'on';
    %     ax.YMinorGrid = 'on';
    %     ax.TickLength = [0.02 0.025];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %IMMX Plots
    %IMMXPlots(parameters,expData,dis,loadn,i,strainRange,figIMMX);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%Plots for the Bristol event
% trousersEventPlots(expData);

%K obtained from the linearization process, when strainRelax is contained
%in the fit range
% save('readInstronData.mat','expData');
% temp = MatPar(1);
% save('MatPar_EPR.mat','-struct','temp');
% temp = MatPar(2);
% save('MatPar_FR.mat','-struct','temp');
% temp = MatPar(3);
% save('MatPar_NatR.mat','-struct','temp');
% temp = MatPar(4);
% save('MatPar_NR.mat','-struct','temp');
% temp = MatPar(5);
% save('MatPar_PR.mat','-struct','temp');
% temp = MatPar(6);
% save('MatPar_SR.mat','-struct','temp');

function [segmentsIndex, segments] = optimizedStrainSegments(stress, strain, error)
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

function trousersEventPlots(expData,figTrousers)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare vectors for ultimate load, toe stiffness and elastic stiffness
loadNames = categorical({'H Tendon';
    'EPR';
    'FR';
    'NatR';
    'NR';
    'PR';
    'SR'});
loadNames = reordercats(loadNames,{'H Tendon';
    'EPR';
    'FR';
    'NatR';
    'NR';
    'PR';
    'SR'});
loadAll = [2172.7;
    expData(1).uLoad;
    expData(2).uLoad;
    expData(3).uLoad;
    expData(4).uLoad;
    expData(5).uLoad;
    expData(6).uLoad];
relaxation = [41;
    14;
    33.2;
    16.8;
    12.5;
    33.3;
    10.3];

figure(figTrousers)
%Switch between double axis plot and single plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
yyaxis left
% bar(loadNames(1),loadAll(1),'EdgeColor',[0.3010,0.7450,0.9330],'FaceColor',[0.3010,0.7450,0.9330]);
bar(loadNames,loadAll,'EdgeColor',[0.3010,0.7450,0.9330],'FaceColor',[0.3010,0.7450,0.9330]);
ylabel('Ultimate Load (N)');

yyaxis right
bar(loadNames(2:end),loadAll(2:end),'EdgeColor',[0.3010,0.7450,0.9330],'FaceColor',[0.3010,0.7450,0.9330]);
ylabel('Ultimate Load (N)');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% bar(loadNames,relaxation,'EdgeColor',[0.3010,0.7450,0.9330],'FaceColor',[0.3010,0.7450,0.9330]);
% ylabel('Relaxation (%)');

grid on;
end

function IMMXPlots(parameters,expData,dis,loadn,i,strainRange,figIMMX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Code for IMMMX paper charts
%Linear Fit for Toe and Elastic Region
elastic = find(dis >= parameters(i).el);
toe = find(dis <= parameters(i).toe);
%     dis = strain;
%     load = stress;
[~,~,stfEl(i),interEl] = linearfit(dis(elastic(1):strainRange(end)),loadn(elastic(1):strainRange(end)));
[~,~,stfToe(i),interToe] = linearfit(dis(toe),loadn(toe));

elFit = stfEl(i)*dis(strainRange) + interEl;
toeFit = stfToe(i)*dis(strainRange) + interToe;
%     parameters(i).toe = stfToe;
%     parameters(i).el = stfEl;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Multiplier calculations IMMX

multFinal(i) = 2172.7./expData(i).uLoad;
newMatAreas(i) = multFinal(i).*parameters(i).area*1e4;

figure(figIMMX)
subplot(2,3,i);
plot(dis(strainRange),loadn(strainRange),'Color',colors(6,:));
hold on
plot(dis(strainRange),toeFit,'--','LineWidth',2,'Color',colors(2,:));
plot(dis(strainRange),elFit,'--','LineWidth',2,'Color',colors(4,:));
hold off
legend('Tensile Strength','Toe Stiffness.','Elastic Stiffness.');
sTitle = parameters(i).name;
title(sTitle);
ylabel('Load (N)');
xlabel('Elongation (mm)');
grid on;
ax = gca; % current axes
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.TickLength = [0.02 0.025];
%     xlim([0 60]);
ylim(parameters(i).yAxisLimits);
end