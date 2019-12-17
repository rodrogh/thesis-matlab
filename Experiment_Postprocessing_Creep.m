%Creep experiment. Type III Specimen
%PR = Polyethylene Rubber
clear all;
close all;
clc;

lineSize = 8;
fontSize = 38;
fileStart = [27 0]; %[Row Column] where data is being held
%Extract data from csv
path = 'RDSO_TypeIII_Creep_1.is_trelax_RawData\\Specimen_RawData_1.csv';
PRData1 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Creep_1.is_trelax_RawData\\Specimen_RawData_2.csv';
PRData2 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Creep_1.is_trelax_RawData\\Specimen_RawData_3.csv';
PRData3 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Creep_1.is_trelax_RawData\\Specimen_RawData_4.csv';
PRData4 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Creep_1.is_trelax_RawData\\Specimen_RawData_5.csv';
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
region = 4051;
gageLength = 50;
width = 19;
thickness = 3;
crossArea = width*1e-3*thickness*1e-3;
stress = 5/crossArea;

plot(PRData1(1:region,1),PRData1(1:region,2)/gageLength,...
     PRData2(1:region,1),PRData2(1:region,2)/gageLength,...
     PRData3(1:region,1),PRData3(1:region,2)/gageLength,...
     PRData4(1:region,1),PRData4(1:region,2)/gageLength,...
     PRData5(1:region,1),PRData5(1:region,2)/gageLength,'LineWidth',lineSize); 
title('Polyethylene Rubber Creep Curve','FontSize',fontSize);
ylabel('Strain (%)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5','FontSize',fontSize);
xlim([0 400]);
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
%Calculation of moving avarage values
%TIME COLUMN IS AVARAGED, IS BETTER TO USE THE NON AVERAGED ONE
k1 = 10; %Moving avarage parameter for elements to include. Center current element
PRDataMean1 = movmean(PRData1,k1,1,'Endpoints','shrink');
PRDataMean2 = movmean(PRData2,k1,1,'Endpoints','shrink');
PRDataMean3 = movmean(PRData3,k1,1,'Endpoints','shrink');
PRDataMean4 = movmean(PRData4,k1,1,'Endpoints','shrink');
PRDataMean5 = movmean(PRData5,k1,1,'Endpoints','shrink');
figure
plot(PRData1(1:region,1), PRDataMean1(1:region,2)/gageLength,...
     PRData2(1:region,1), PRDataMean2(1:region,2)/gageLength,...
     PRData3(1:region,1), PRDataMean3(1:region,2)/gageLength,...
     PRData4(1:region,1), PRDataMean4(1:region,2)/gageLength,...
     PRData5(1:region,1), PRDataMean5(1:region,2)/gageLength,'LineWidth',lineSize);
%plotTitle=sprintf('Polyethylene Rubber Time vs Stress - %d points Moving Average',k1);
plotTitle='Polyethylene Rubber Creep';
title(plotTitle,'FontSize',fontSize);
ylabel('Strain (%)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5','FontSize',fontSize);
%axis([0 region/10 0 inf])
xlim([0 400]);
grid on;
ax = gca; % current axes
ax.FontSize = fontSize;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.TickLength = [0.02 0.025];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Unifying all specimen's data into one
PRDataDisc = zeros(region+1,6);
%The range is from 2:end to have an initial zero as the first element
PRDataDisc(2:region+1,:) = PRData1(1:region,:) + PRData2(1:region,:) +...
                           PRData3(1:region,:) + PRData4(1:region,:) +...
                           PRData5(1:region,:);
PRDataDisc = PRDataDisc/5;

figure
plot(PRDataDisc(:,1),PRDataDisc(:,2)/gageLength);
plotTitle = sprintf('Polyethylene Rubber Time vs Strain - Unified chart');
title(plotTitle,'FontSize',fontSize);
ylabel('Strain (%)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
axis([0 region/10 0 inf])
ax = gca; % current axes
ax.FontSize = fontSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stiffness: diff() to obtain approximate derivate
%Approximate derivative with diff()
%2nd parameter = order
%3rd parameter = dimension along to perform differentiation
k2 = 5;
PRDataDiscDer = zeros(region+1,6);
PRDataDiscDer(2:end,:) = diff(PRDataDisc/gageLength,1,1);

%Stress Relaxation rate = d(load)/d(time)
PRDataStressRelaxation = PRDataDiscDer(:,2)./PRDataDiscDer(:,1);
PRDataStressRelaxation(1,1) = 0;

%Plot Stress Relaxation
%The experiment deformation rate is 1mm/sec
%Samples are available for each 100 ms. Therefore first 200 data points are
%plotted
figure
plot(PRDataDisc(1:region,1),PRDataStressRelaxation(1:region,1),'r--');
plotTitle = sprintf('Polyethylene Rubber Time vs Strain Rate - Unified chart');
title(plotTitle,'FontSize',fontSize);
ylabel('Strain Rate (%/s)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
axis([0 region/10 0 inf])
ax = gca; % current axes
ax.FontSize = fontSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the creep compliance and the relaxation time
Ccrp = (PRDataDisc(:,2)/gageLength)/stress;
figure
semilogx(PRDataDisc(:,1),Ccrp(:),'LineWidth',lineSize);
%plotTitle='Polyethylene Rubber Creep Compliance';
plotTitle='Ensayo de Fluencia';
title(plotTitle,'FontSize',fontSize);
ylabel('Modulo de Fluencia','FontSize',fontSize);
xlabel('Tiempo (segundos)','FontSize',fontSize);
grid on;
ax = gca; % current axes
ax.FontSize = fontSize;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.TickLength = [0.02 0.025];
figure
semilogx(PRDataDisc(1:end-2,1),diff(diff(Ccrp)))





