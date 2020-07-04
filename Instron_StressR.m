%Postprocessing of Instron Stress Relaxation Experiment data
%The first important feature to obtain is the Yield Strength of the
%material in order to stablish a linear region where the material should be
%when performing time-dependant experiments

clear all;
close all;
clc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The data is distributed in six columns as follows:
%1 - Time
%2 - Extension
%3 - Load
%4 - Tensile Strain
%5 - Tensile Stress
%6 - Corrected position
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Font parameters for plotting
lineSize = 2;
fontSize = 18;
font = 'Gill Sans MT';
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lineSize);   % set the default line width to lw
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

figStress = figure('Name','Stress Relaxation');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameter for data extraction. Acces each element with filePaths(n).paths
lo = 33; %initial length
parameters = struct('name',...
    {'Ethylene Polypropylene Rubber';
    'Fluorocarbon Rubber';
    'Natural Rubber';
    'Nitrile Rubber';
    'Polyethylene Rubber';
    'Silicone Rubber';},...
    'fieldName',...
    {'EPR';
    'FR';
    'NatR';
    'NR';
    'PR';
    'SR';},...
    'paths',...
    {'ASTM D412\Stress Relaxation\RDSO EPR SRelax 500.is_trelax_RawData\Specimen_RawData_';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 500.is_trelax_RawData\Specimen_RawData_';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 500.is_trelax_RawData\Specimen_RawData_';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 500.is_trelax_RawData\Specimen_RawData_';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 500.is_trelax_RawData\Specimen_RawData_';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 500.is_trelax_RawData\Specimen_RawData_'},...
    'area',...
    {(6e-3*1.5e-3);    (6e-3*1.5e-3);    (6e-3*1.5e-3);    (6e-3*1.5e-3);  (6e-3*6e-3);  (6e-3*1.5e-3)},...
    'branches',...
    {9; 8; 10; 8; 8; 10});
%Original number of branches
% 'branches',...
%     {4; 8; 4; 6; 6; 6});
%Number of branches found with the alternative method (High) - Eng viscoelasticity
% 'branches',...
%     {11; 11; 11; 13; 13; 10});
%Number of branches found with the alternative method (Low) - Eng viscoelasticity
% 'branches',...
%     {9; 8; 10; 8; 8; 10});
expData = struct('data',{},'time',{},'strain',{},'stress',{},'YMtoe',{},'YM',{},...
    'uStrain',{},'uStress',{},'uDis',{},'uLoad',{},'load',{},'dis',{},...
    'stiffToe',{},'stiffEl',{});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Call to ReadInstronTable is not required now. Files have been created
%compiling the processed data from all tests.
%SRelSetProc.mat - Averaged (per experimental setup) data structure
%SRelSetProcInd.mat - Individually processed tests
%Method:
%0 for normal processing behaviour
%1 for neural network processing behaviour
method = 0;
if method == 0
    load('../Experimental Work/TrainingSetProcSRelSetProcStressV4.mat');
else
    load('../Experimental Work/SRelSetProcInd.mat'); %This approach is probably unnecesary
end
%Core variables initialization
%Structure to save fitted data for the all the different
%experimental setups
fitParameters = PmatData;
%Defining number of materials and number of experimental setups
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR';...
%RBYellow, RBRed, RBBlue, RBGreen, RBBlack, RBOrange}
fields.N = fieldnames(PmatData);
M = length( fieldnames( PmatData.(fields.N{1}) ) );
%Exp setups NAMES: {'L5';'L15';'L180'}
fields.M = fieldnames( PmatData.(fields.N{1}) );
%Data variables NAME {'paths';'pload';'rload';'rdis';'rtime'}
P = length( fieldnames( PmatData.(fields.N{1}).(fields.M{1}) ) );
fields.P = fieldnames( PmatData.(fields.N{1}).(fields.M{1}) );

k = zeros(N,6+1); %6 is the maximum number of branches
eta = zeros(N,6);
thao = zeros(N,6);
Eo = zeros(N,1);
%Create gof variable and fill it with nan to avoid zero values in the
%min calculation
Branches = 10;
gofSR = nan(N,Branches);

for i=1:N
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    M = length(fields.M);
    for j=1:M  %{'L5';'L15';'L180'}
        %Initialize structure to save fitted parameters
        fitParameters.( fields.N{i} ).( fields.M{j} ) = [];
        %Create more ergonomic variables for each parameter
        load = PmatData.( fields.N{i} ).( fields.M{j} ). sstress;
        dis = PmatData.( fields.N{i} ).( fields.M{j} ). pstrain;
        time = PmatData.( fields.N{i} ).( fields.M{j} ). ptime;
        %Experimental setup: Initial deformation
        Eo(i) = dis(end);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %SLS model fit for sress relaxation:
        [kSLS(i,:), etaSLS(i,:), thaoSLS(i,:), EoSLS(i,:)] = SLSStressRel(load,Eo(i),time);
%         [kSLS(i,:), etaSLS(i,:), thaoSLS(i,:), EoSLS(i,:)] = optimizingSLS(load,Eo(i),time);
        
        %Applying formula deducted from exponential decaying curve
        loadSLS = EoSLS(i)*(kSLS(i,1) + kSLS(i,2)*(1./exp(time./thaoSLS(i))));
%         loadSLS2 = EoSLS2(i)*(kSLS2(i,1) + kSLS2(i,2)*(1./exp(time./thaoSLS2(i))));
        gofSLS(i) = sqrt( mean( (load-loadSLS).^2) );
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Wiechert model fit with 2 branches
        %     [kW(i,:), etaW(i,:), thaoW(i,:), EoW(i,:)] = wiechertExp2StressRel(stress,Eo(i),time);
        %     [kW(i,:), etaW(i,:), thaoW(i,:), EoW(i,:)] = wiechertStressRel(stress,Eo(i),time,2,parameters(i));
        %     [kW(i,:), etaW(i,:), thaoW(i,:), EoW(i,:)] = optimizingWiechert2(load,Eo(i),time,2,parameters(i));
        %
        %     stressW2 = EoW(i)*(kW(i,1) + kW(i,2)*(1./exp(time./thaoW(i,1))) + kW(i,3)*(1./exp(time./thaoW(i,2))));
        %     gofW2(i) = 100*(1/mean(load))*sqrt(mean(power(load-stressW2,2)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Wiechert model fit with optimized  N branches
        %Define k as a cell matrix of NxM dimension which holds a structure
        %araay with all the branches tested
        tempCell = cell(N,M);
        fitOpt = struct([]); %Fit Optimization structure
        stressRel = zeros(Branches,length(time));
        gofSR = zeros(1,Branches);
        %This for loop was designed to test the effect of the quantity of
        %branches iteratively. If the fit of a particular single value of
        %branches is desired, the for loop can be bipassed making it "for
        %branches=1".        
        for branches=1:Branches
                    [fitOpt(branches).k,fitOpt(branches).eta,fitOpt(branches).thao,Eo(i)] = ...
                        wiechertStressRel(load,Eo(i),time,branches,1);
                    %This variable changes in every iteration
                    tempk = fitOpt(branches).k; %1xbranches+1 vector
                    tempThao = fitOpt(branches).thao; %1xbranches vector
                    %The stress(t) = k0*Eo + SUMoF_j_To_branches kj*exp(-t/thaoj)
                    for test=1:length(time)
                       stressRel(branches,test) =  Eo(i)*tempk(1) + ... %Equilibrium Spring Stiffness
                           Eo(i)*sum( ...
                           tempk(2:end)./exp(time(test)./tempThao)... %This operation results in a 1xbranches vector
                           ); 
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %Compare mathematical model with experimental data and exponential fit
                    gofSR(branches) = sqrt( mean( (load-stressRel(branches,:)).^2) );
        end
        % Calculate optimal number of branches for wiechert model        
        [minSR,indexSR] = min(gofSR,[],2);
        bestStress = stressRel(indexSR,:);
        bestBranches = indexSR;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%         colors = get(gca,'ColorOrder');
%         %         subplot(2,3,i);
%         plot(time,load,'Color',colors(6,:),...
%             'DisplayName','Experimental Data');
%         hold on
%         plot(time,bestStress,'--','Color',colors(7,:),...
%             'DisplayName','Wiechert Fit');
%         hold off
%         legend('Location','northeast','Interpreter','latex');
%         sTitle = strcat("Stress Relaxation Fit - ",fields.N{i}," - ", fields.M{j}," - ",num2str(indexSR));
%         title(sTitle,'Interpreter','latex');
%         xlabel('Time (minutes)');
%         ylabel('Stress (MPa)');
%         %grid on;
%         ax = gca; % current axes
%         ax.Box = 'off';
%         ax.XMinorTick = 'on';
%         ax.YMinorTick = 'on';
%         ax.XGrid = 'on';
%         ax.YGrid = 'on';
%         ax.TickLength = [0.02 0.025];
% %         ax.XLim = [0 180];
% %         ax.YLim = [0 1.1];
        
        %Structure to save fitted data for the all the different
        %experimental setups
        fitParameters.( fields.N{i} ).( fields.M{j} ).SLS.k = kSLS(i,:);
        fitParameters.( fields.N{i} ).( fields.M{j} ).SLS.thao = thaoSLS(i);
        fitParameters.( fields.N{i} ).( fields.M{j} ).SLS.eta = etaSLS(i);
        fitParameters.( fields.N{i} ).( fields.M{j} ).SLS.Eo = Eo(i);
        fitParameters.( fields.N{i} ).( fields.M{j} ).W.k = fitOpt(indexSR).k;
        fitParameters.( fields.N{i} ).( fields.M{j} ).W.thao = fitOpt(indexSR).thao;
        fitParameters.( fields.N{i} ).( fields.M{j} ).W.eta = fitOpt(indexSR).eta;
        fitParameters.( fields.N{i} ).( fields.M{j} ).W.bestBranches = indexSR;
        fitParameters.( fields.N{i} ).( fields.M{j} ).W.Eo = Eo(i);
    end
end

% Save Wiechert model fit
if method == 0
    save('../Experimental Work/SRelFitParametersV5.mat','fitParameters');
else
    save('../Experimental Work/SRelFitParametersIndV5.mat','fitParameters');
end
% save('readInstronDataSR.mat','expData');
% save('readInstronRawDataSR.mat','expData');

function [k,eta,thao,Eo] = wiechertStressRel(stress,Eo,time,N,parameters)
%Find the coefficients ki and etai from N Maxwell branches during a stress
%relaxation test where the strain is constant during the test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Algorithm to spread the point of interest between the decades 10^1 - 10^3
%of the time vector. Need to find the point in time where the first decade
%starts and where the third decades end. Between this range we want to
%logspace our points of interest
%Convert time to logarithmic scale
timeLog = log10(time);
%Create a linearly spaced sampling vector. Discarding first and last
%datapoints because they can't be used as time constants
%Start from timeLog(2), instead of timeLog(1)
%End at timeLog(end-1), instead of timeLog(end)
timeLog_LinScale = linspace(timeLog(2),timeLog(end-1),N);
indexOfTimeCollocation = zeros(1,N);
for i=1:N
    indexes = find( timeLog >= timeLog_LinScale(i));
    indexOfTimeCollocation(i) = indexes(1);
end
%Extract N time constants
thao = time(indexOfTimeCollocation); %1xN vector
stressCollocation = stress(indexOfTimeCollocation);
%Obtain the Relaxation Modulus E(t) = sigma(t)/Eo; inclluding t=0.
%This will result in a 1x(N+1) vector
E = [stress(1)/Eo stressCollocation/Eo]'; %Nx1 vector
%% Time Collocation approach by C. Machiraju
%At t=0 the sum of all k values is equal to E(t=0)
%This equation is useful for the calculation and is in addition to the N
%time constants desired to collocate.
%Therefore the equation system to create has N+1 equations, as follows:
%Define system of equations relating each time constant to a known
%relaxation modulus (E)
A = ones(N+1,N+1);
%The equation is solved by doing A*k = E; -> k=inv(A)*E;
%Start from [2,2] to [N+1,N+1]
for i=2:N+1
    for j=2:N+1
        A(i,j) = 1./exp(thao(i-1)/thao(j-1));
    end
end
%Solve
k = A\E;
%Obtain viscous coefficients
eta = thao.*(k(2:end)');
%Obtain the transpose of all variables to match calling function
k = k';
% figure
% semilogx(time_i,stress_i*1e-6,'o');
% legend('Collocation Log');
%         sTitle = strcat(parameters.name,' - Time Logarithmic Collocation');
%         title(sTitle);
%         xlabel('Time (seconds)');
%         ylabel('Stress (MPa)');
%         %grid on;
%         ax = gca; % current axes
%         ax.XMinorTick = 'on';
%         ax.YMinorTick = 'on';
%         ax.XMinorGrid = 'on';
%         ax.YMinorGrid = 'on';
%         ax.TickLength = [0.02 0.025];
end

function [k,eta,thao,Eo] = SLSStressRel(stress,Eo,time)
%Equilibrium spring's stiffness
ke = stress(end)/Eo;
km = max(stress)/Eo - ke;
stress_at_thao = Eo*(ke + km/exp(1));
thaoIndex = find(stress > stress_at_thao);
thao = time(thaoIndex(end));
eta = thao*km;
k = [ke km];
end

function [k, eta, thao, Eo] = wiechertExp2StressRel(stress,Eo,time)
%Fitting the 2 exponential eq a*exp(b*x) + c*exp(d*x) requires an extra
%step since this Eq does not account for the stress from ke. The actual
%eq for a 2-term wiechert is S = (ke + a*exp(b*x) + c*exp(d*x))*Eo. Before
%fitting I need to find (S/Eo) - ke = a*exp(b*x) + c*exp(d*x)
ke = stress(end)/Eo;
stressFit = stress - stress(end);
% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time, stressFit );
% Set up fittype and options.
ft = fittype( 'exp2' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Upper = [Inf 0 Inf 0];
opts.Lower = [-Inf -time(end) -Inf -time(end)];
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
% % Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'stress vs. time', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel time
% ylabel stress
% grid on
%Following a*exp(b*x) + c*exp(d*x)
%Extract thaos and calculate K's using equation system
coeff = coeffvalues(fitresult);
thao(1) = -1/coeff(2);
thao(2) = -1/coeff(4);
N=2;
for n=1:N
    range = find(time < thao(n));
    if isempty(range)
        points(n) = 1;
        thao(n) = time(1);
    else
        points(n) = range(end);
    end
end
stress_i = stress(points);
%Equation system only solve the variables of the branches
A2 = zeros(N,N);
B2 = zeros(N,1);
for i=1:N
    B2(i,1) = stress_i(i)/Eo - ke;
    for j=1:N
        A2(i,j) = 1./exp(thao(i)/thao(j));
    end
end
%Choose which method to use
%Solve equation system by finding k = inv(A)*B
k = A2\B2;
k(1) = stress(1)/Eo - ke - sum(k(2:end));
k = [ke ; k];
temp = k';
clear k
k = temp;
eta = thao.*k(2:end);
end

function [k,eta,thao,Eo] = optimizingWiechert2(stress,Eo,time,N,parameters)
%Everything is tailored for 2 branches
points = ones(1,4);
points(1) = 1;
points(end) = length(time);
a = 1;
b = 1;
for thao1 = 10:100:time(end)
    % for thao1 = 10 + (5-1)*100
    range = find(time < thao1);
    points(2) = range(end);
    for thao2 = 50:100:time(end)
        %     for thao2 = 50 + (42-1)*100
        range = find(time < thao2);
        points(3) = range(end);
        time_i = time(points);
        stress_i = stress(points);
        thao = time_i(2:3);
        matrixthao(a,b,:) = thao; %Holds all iterations
        %Alternative method - Enginnering viscoelasticity
        %Equation system only solve the variables of the branches
        A2 = zeros(N,N);
        B2 = zeros(N,1);
        ke = stress(end)/Eo;
        for i=1:N
            B2(i,1) = stress_i(i+1)/Eo - ke;
            for j=1:N
                A2(i,j) = 1./exp(thao(i)/thao(j));
            end
        end
        k = A2\B2;
        k(1) = stress(1)/Eo - ke - sum(k(2:end));
        k = [ke ; k];
        matrixk(a,b,:) = k; %Holds all iterations
        %%%%%%%%%%%
        %Obtain viscous coefficients
        eta = thao.*k(2:end);
        matrixeta(a,b,:) = eta;%Holds all iterations
        %obtain the transpose of all variables to match calling function
        temp = k';
        clear k
        k = temp;
        temp = eta';
        clear eta
        eta = temp;
        temp = thao';
        clear thao
        thao = temp;
        %Calculate stress response
        stressW2 = Eo*(k(1) + k(2)*(1./exp(time./thao(1))) + k(3)*(1./exp(time./thao(2))));
        
        %Compare mathematical model with experimental data and exponential fit
        gof(a,b) = mean(sqrt(power(stress-stressW2,2)));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                 figure
        %                 colors = get(gca,'ColorOrder');
        %                 plot(time/60,stress*1e-6,'Color',colors(6,:))
        %                 hold on
        %                 plot(time/60,stressW2*1e-6);
        %                 hold off
        %                 legend('Data','WiechertLog2 Fit');
        %                 sTitle = strcat(parameters.name,'');
        %         %         sTitle = strcat(parameters(i).name,' - Stress Relaxation');
        %                 title(sTitle);
        %                 xlabel('Time (minutes)');
        %                 ylabel('Stress (MPa)');
        %                 %grid on;
        %                 ax = gca; % current axes
        %                 ax.XMinorTick = 'on';
        %                 ax.YMinorTick = 'on';
        %                 ax.XMinorGrid = 'on';
        %                 ax.YMinorGrid = 'on';
        %                 ax.TickLength = [0.02 0.025];
        %                 ax.XLim = [0 180];
        %                 %ax.YLim = [0 1.1];
        b = b + 1;
    end
    b = 1;
    a = a + 1;
end
%optimal number of branches = 5th iteration in thao1 and 42th iteration in
%thao2. The index for these values are: 202 and 2072
[mincol,indexcol] = min(gof,[],2);
[minrow,indexrow] = min(mincol);

temp = matrixthao(indexrow,indexcol(indexrow),:);
thao = temp(1,:)';
temp = matrixeta(indexrow,indexcol(indexrow),:);
eta = temp(1,:)';
temp = matrixk(indexrow,indexcol(indexrow),:);
k = temp(1,:)';
end

function [k,eta,thao,Eo] = optimizingSLS(load,Eo,time)
%This function does not work as it is
%Equilibrium spring's stiffness
ke = load(end)/Eo;
km = max(load)/Eo - ke;
load_at_thao = Eo*(ke + km/exp(1));
thaoSpace = linspace( time(1) , time(end) , 1000);

for i = 1:length(thaoSpace)
    loadFit = Eo*(ke + km./exp(time./thaoSpace(i)));
    gof(i) = immse( load , loadFit) ;
end

[best_gof,i_m] = min(gof);
thao = thaoSpace(i_m);
eta = thao*km;
k = [ke km];
end