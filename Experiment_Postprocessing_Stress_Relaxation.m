%Stress relaxation experiment. Type III Specimen
%PR = Polyethylene Rubber
clear all;
close all;
clc;

lineSize = 6;
fontSize = 38;
fileStart = [15 0]; %[Row Column] where data is being held
%Extract data from csv
path = 'RDSO_TypeIII_Stress_Relaxation_1.is_trelax_RawData\\Specimen_RawData_1.csv';
PRData1 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Stress_Relaxation_1.is_trelax_RawData\\Specimen_RawData_2.csv';
PRData2 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Stress_Relaxation_1.is_trelax_RawData\\Specimen_RawData_3.csv';
PRData3 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Stress_Relaxation_1.is_trelax_RawData\\Specimen_RawData_4.csv';
PRData4 = csvread(path,fileStart(1),fileStart(2));
path = 'RDSO_TypeIII_Stress_Relaxation_1.is_trelax_RawData\\Specimen_RawData_5.csv';
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
region = 4033;
gageLength = 50;
width = 19;
thickness = 3;
crossArea = width*1e-3*thickness*1e-3;
strain = 0.06;
% subplot(2,1,1); plot(PRData1(1:region,1),PRData1(1:region,3)/crossArea,...
%      PRData2(1:region,1),PRData2(1:region,3)/crossArea,...
%      PRData3(1:region,1),PRData3(1:region,3)/crossArea,...
%      PRData4(1:region,1),PRData4(1:region,3)/crossArea,...
%      PRData5(1:region,1),PRData5(1:region,3)/crossArea); 
% title('Polyethylene Rubber Time vs Stress - Raw Data','FontSize',fontSize);
% ylabel('Stress (N/m^2)','FontSize',fontSize);
% xlabel('TIme (seconds)','FontSize',fontSize);
% legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5','FontSize',fontSize);
% axis([0 region/10 0 inf])
% ax = gca; % current axes
% ax.FontSize = fontSize;
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
%subplot(2,1,2); 
plot(PRData1(1:region,1), PRDataMean1(1:region,3)/crossArea,...
     PRData2(1:region,1), PRDataMean2(1:region,3)/crossArea,...
     PRData3(1:region,1), PRDataMean3(1:region,3)/crossArea,...
     PRData4(1:region,1), PRDataMean4(1:region,3)/crossArea,...
     PRData5(1:region,1), PRDataMean5(1:region,3)/crossArea,'LineWidth',lineSize);
%plotTitle=sprintf('Polyethylene Rubber Time vs Stress - %d points Moving Average',k1);
plotTitle='Polyethylene Rubber Stress Relaxation';
title(plotTitle,'FontSize',fontSize);
ylabel('Stress (N/m^2)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
legend('Specimen 1','Specimen 2','Specimen 3','Specimen 4','Specimen 5','FontSize',fontSize);
%axis([0 region/10 0 inf])
xlim([0 400]);
ylim([0 6.5e4]);
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
PRDataDisc = PRDataMean1 + PRDataMean2 +...
                           PRDataMean3 + PRDataMean4 +...
                           PRDataMean5;
PRDataDisc = PRDataDisc/5;

figure
plot(PRDataDisc(:,1),PRDataDisc(:,3)/crossArea,'LineWidth',lineSize);
plotTitle = sprintf('Polyethylene Rubber Stress Relaxation Unified');
title(plotTitle,'FontSize',fontSize);
ylabel('Stress (N/m^2)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
xlim([0 400]);
ylim([0 5.5e4]);
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
%Stiffness: diff() to obtain approximate derivate
%Approximate derivative with diff()
%2nd parameter = order
%3rd parameter = dimension along to perform differentiation
k2 = 5;
PRDataDiscDer = zeros(4033,6);
PRDataDiscDer(2:end,:) = diff(PRDataDisc/crossArea,1,1);
PRDataDiscDerMean = movmean(PRDataDiscDer,k2,1);

%Stress Relaxation rate = d(load)/d(time)
PRDataStressRelaxation = PRDataDiscDer(:,3)./PRDataDiscDer(:,1);
PRDataStressRelaxation(1,1) = 0;
PRDataStressRelaxationMean = PRDataDiscDerMean(:,3)./PRDataDiscDerMean(:,1);
% PRDataStiffnessToe = sum(PRDataStressRelaxation(toeRegion))/length(toeRegion);
% PRDataStiffnessElastic = sum(PRDataStressRelaxation(elasticRegion))/length(elasticRegion);
%Plot Stress Relaxation
%The experiment deformation rate is 1mm/sec
%Samples are available for each 100 ms. Therefore first 200 data points are
%plotted
figure
subplot(2,1,1); 
plot(PRDataDisc(1:region,1),PRDataStressRelaxation(1:region,1),'r--');
plotTitle = sprintf('Polyethylene Rubber Time vs Stress Relaxation Rate - Unified chart');
title(plotTitle,'FontSize',fontSize);
ylabel('Stress Relaxation Rate (Pa/s)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
axis([0 region/10 0 inf])
ax = gca; % current axes
ax.FontSize = fontSize;

subplot(2,1,2); 
plot(PRDataDisc(1:region,1),PRDataStressRelaxationMean(1:region,1),'b');
plotTitle = sprintf('Polyethylene Rubber Time vs Stress Relaxation Rate - %d points Moving Average',k2);
title(plotTitle,'FontSize',fontSize);
ylabel('Stress Relaxation Rate (Pa/s)','FontSize',fontSize);
xlabel('Time (seconds)','FontSize',fontSize);
axis([0 region/10 0 inf])
ax = gca; % current axes
ax.FontSize = fontSize;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the relaxation modulus and the relaxation time
Erel = (PRDataDisc(:,3)/crossArea)/strain;
figure
semilogx(PRDataDisc(:,1),Erel(:),'LineWidth',lineSize);
%plotTitle='Polyethylene Rubber Stress Rel. Modulus';
plotTitle='Ensayo de Relajacion de Tensiones';
title(plotTitle,'FontSize',fontSize);
ylabel('Modulo de Relajacion (N/m^2)','FontSize',fontSize);
xlabel('Tiempo (segundos)','FontSize',fontSize);
grid on;
ax = gca; % current axes
ax.FontSize = fontSize;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.XMinorGrid = 'on';
ax.YMinorGrid = 'on';
ax.TickLength = [0.02 0.025];


