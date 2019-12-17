%Hysteresis experiment. Type III Specimen
%PR = Polyethylene Rubber
clear all;
close all;
clc;

lineSize = 8;
fontSize = 38;
fileStart = [13 0]; %[Row Column] where data is being held
%Extract data from csv
path = 'RDSO_TypeIII_Hysteresis_1mmsec.is_tcyclic_RawData\\Specimen_RawData_1.csv';
PRData1 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Hysteresis_1mmsec.is_tcyclic_RawData\\Specimen_RawData_2.csv';
PRData2 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Hysteresis_1mmsec.is_tcyclic_RawData\\Specimen_RawData_3.csv';
PRData3 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Hysteresis_1mmsec.is_tcyclic_RawData\\Specimen_RawData_4.csv';
PRData4 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Hysteresis_1mmsec.is_tcyclic_RawData\\Specimen_RawData_5.csv';
PRData5 = csvread(path,fileStart(1),fileStart(2));
%Array data per column:
%C1 = Time
%C2 = Extenstion
%C3 = Load
%C4 = Tensile Extension
%C5 = Tensile Stress
%C6 = Tensile Strain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot of raw data
figure 
% subplot(2,1,1); plot(PRData1(:,2),PRData1(:,3),...
%      PRData2(:,2),PRData2(:,3),...
%      PRData3(:,2),PRData3(:,3),...
%      PRData4(:,2),PRData4(:,3),...
%      PRData5(:,2),PRData5(:,3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Region has to be smaller than the smallest PRData matrix
region = 102;
subplot(2,1,1); plot(PRData1(1:region,2),PRData1(1:region,3),...
     PRData2(1:region,2),PRData2(1:region,3),...
     PRData3(1:region,2),PRData3(1:region,3),...
     PRData4(1:region,2),PRData4(1:region,3),...
     PRData5(1:region,2),PRData5(1:region,3)); 
title('Polyethylene Rubber Extension vs Load - Raw Data');
ylabel('Load (N)');
xlabel('Extension (mm)');
legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5');

%Calculation of moving avarage values
%TIME COLUMN IS AVARAGED, IS BETTER TO USE THE NON AVERAGED ONE
k1 = 10; %Moving avarage parameter for elements to include. Center current element
PRDataMean1 = movmean(PRData1,k1,1,'Endpoints','shrink');
PRDataMean2 = movmean(PRData2,k1,1,'Endpoints','shrink');
PRDataMean3 = movmean(PRData3,k1,1,'Endpoints','shrink');
PRDataMean4 = movmean(PRData4,k1,1,'Endpoints','shrink');
PRDataMean5 = movmean(PRData5,k1,1,'Endpoints','shrink');
%Plot of moving avarage of Load vs Extension
% subplot(2,1,2); plot(PRDataMean1(:,2), PRDataMean1(:,3),...
%      PRDataMean2(:,2), PRDataMean2(:,3),...
%      PRDataMean3(:,2), PRDataMean3(:,3),...
%      PRDataMean4(:,2), PRDataMean4(:,3),...
%      PRDataMean5(:,2), PRDataMean5(:,3));
subplot(2,1,2); plot(PRData1(1:region,2), PRDataMean1(1:region,3),...
     PRData2(1:region,2), PRDataMean2(1:region,3),...
     PRData3(1:region,2), PRDataMean3(1:region,3),...
     PRData4(1:region,2), PRDataMean4(1:region,3),...
     PRData5(1:region,2), PRDataMean5(1:region,3));
plotTitle=sprintf('Polyethylene Rubber Extension vs Load - %d points Moving Average',k1);
title(plotTitle);
ylabel('Load (N)');
xlabel('Extension (mm)');
legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unifying all specimen's data into one
%Discretizing
PRDataDisc = zeros(region+1,6);
%The range is from 2:end to have an initial zero as the first element
PRDataDisc(2:region+1,:) = PRDataMean1(1:region,:) + PRDataMean2(1:region,:) +...
                           PRDataMean3(1:region,:) + PRDataMean4(1:region,:) +...
                           PRDataMean5(1:region,:);
PRDataDisc = PRDataDisc/5;

figure
plot(PRDataDisc(:,2),PRDataDisc(:,3),'LineWidth',lineSize);
plotTitle = sprintf('Polyethylene Rubber Extension vs Load - Unified Data');
%plotTitle = 'Ensayo de Histeresis';
title(plotTitle);
ylabel('Load (N)');
xlabel('Strain (mm)');
grid on
ax = gca; % current axes
ax.FontSize = fontSize;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.TickLength = [0.02 0.025];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stiffness: diff() to obtain approximate derivate
%Approximate derivative with diff()
%2nd parameter = order
%3rd parameter = dimension along to perform differentiation
k2 = 5;
PRDataDiscDer = zeros(region+1,6);
PRDataDiscDer(2:end,:) = diff(PRDataDisc,1,1);
PRDataDiscDerMean = movmean(PRDataDiscDer,k2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stress-Strain calculation
%Specimen dimensions (mm)
gageLength = 50;
width = 19;
thickness = 3;
crossArea = width*1e-3*thickness*1e-3; %meters square
