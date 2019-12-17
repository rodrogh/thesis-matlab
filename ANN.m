%MATLAB neural networks tool. Shallow networks for prediction
%Define data holder
clear all;
close all;
clc;
dataStruct = struct('load',{},'pload',{},'dis',{},'time',{});
%Defining paths
% paths = {'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_8';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_9';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_10';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_11';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_12';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_13';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 500.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 500.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\SR\RDSO SR Tensile 500.is_tens_RawData\Specimen_RawData_3'};
paths = {'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 250.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile Multi.is_tens_RawData\Specimen_RawData_1'};
% paths = {'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_8';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_9';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_8';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile Multi.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\FR\RDSO FR Tensile Multi.is_tens_RawData\Specimen_RawData_2'};
% paths = {'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_8';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_9';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile Multi.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NatR\RDSO NatR Tensile Multi.is_tens_RawData\Specimen_RawData_2'};
% paths = {'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile Multi.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\NR\RDSO NR Tensile Multi.is_tens_RawData\Specimen_RawData_2'};
% paths = {'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_8';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_9';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_10';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_11';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_12';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_2';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_3';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_4';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_5';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_6';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_7';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile Multi.is_tens_RawData\Specimen_RawData_1';
%     'ASTM D412\Tensile\PR\RDSO PR6 Tensile Multi.is_tens_RawData\Specimen_RawData_2'};
%Define batch to read
N = length(paths);
iRow = 2; %initial row
iCol = 0; %initial column
%Variable to hold the length of each file
a = zeros(1,N);
%Read data from data tables
%Specimens dimensions to convert load-dis to stress-strain
lo = 33; %Initial length
crossArea = 1.5e-3*1.5e-3;
for i=1:N
    p = strcat(paths{i},'.csv');
    temp = csvread(p,iRow,iCol); %Read rows and columns from file
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %The data is distributed in six columns as follows:
    %1 - Time
    %2 - Extension
    %3 - Load
    %4 - Tensile Strain
    %5 - Tensile Stress
    %6 - Corrected position
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Save desired raw data
    %Find maximum load, before failure
    [maxLoad,maxIndex] = max(temp(:,3));
    dataStruct(i).load = temp(1:maxIndex,3);
    dataStruct(i).dis = temp(1:maxIndex,2);
    dataStruct(i).time = temp(1:maxIndex,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Apply the same preprocessing algorithm to the raw data before training
    %the neural network
%     dataStruct(i).pload = smooth(temp(1:maxIndex,3),0.1,'rloess'); %smooth also filters the NaN values
%     if dataStruct(i).pload(1) > 0
%         dataStruct(i).pload = dataStruct(i).pload - dataStruct(i).pload(1);
%     end    
    [row,~] = size(dataStruct(i).load);
    a(i) = row; %Save the length of each file (rows)
end
totalRow = sum(a);
%The network requires two inputs: the strain and the strain rate
iVar = 2;
% iVar = 1;
dataInput = zeros(totalRow,iVar); %Create input matrix
dataOutput = zeros(totalRow,1);
% dataInputP = zeros(totalRow,iVar); %Create input matrix for preprocess data
% dataOutputP = zeros(totalRow,1);
%Fill input and output matrices
s = 1;
e = a(1);
for i=1:N
    
    dataInput(s+1:e,2) = diff((dataStruct(i).dis))./diff(dataStruct(i).time);
    dataInput(s:e,1) = dataStruct(i).dis;
    %Load corrections for unbalanced load cell. The initial load value is
    %substracted from the whole dataset to achieve 0 load as the initial
    %value and get rid of the unbalanced load offset
    dataOutput(s:e,1) = dataStruct(i).load - dataStruct(i).load(1);
    
%     dataInputP(s:e-1,2) = dataInput(s:e-1,2);
%     dataInputP(s:e,1) = dataInput(s:e,1);
%     dataOutputP(s:e,1) = dataStruct(i).pload;
    
    s = e + 1;
    if i < N
        e = s + a(i+1) - 1;
    end
end
%Save the data which is relevant to the application deformation range of 0
%to 25 mm, i.e. initial sections of each experiment
%Deformation expected in real applications (mm)
% disApp = 25;
% curveSection = find(dataInputP(:,1) < disApp);
% %Raw data
% dataInputApp = dataInput(curveSection,:);
% dataOutputApp = dataOutput(curveSection);
% %Pre-processed data
% dataInputAppP = dataInputP(curveSection,:);
% dataOutputAppP = dataOutputP(curveSection);
%Shift inputs and output to account for previous values
%Created shifted values of both inputs: the displacement and displacement
%rate
temp = dataInput;
% temp = dataInputP;
s1 = circshift(temp,[1 0]); %Shift both dimensions by 1
s1(1,:) = 0;
s2 = circshift(temp,[2 0]); %Shift both dimensions by 2
s2(1:2,:) = 0;
%Save shifted values in order, first displacement and its t,t-1,t-2
inputShift = [temp(:,1),s1(:,1),s2(:,1),...
                temp(:,2),s1(:,2),s2(:,2)];
temp = dataOutput;
% temp = dataOutputP;
s1 = circshift(temp,1);
s1(1) = 0;
s2 = circshift(temp,2);
s2(1:2) = 0;
outputShift = [temp,s1,s2]; %Save shifted values t,t-1,t-2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Saving the formatted dataset
save('NNDataset_Raw_EPR.mat','dataInput','dataOutput',...
   'inputShift','outputShift');
% save('NNDataset_P_PR6.mat','dataInputP','dataOutputP','dataInputApp',...
%     'dataOutputApp','inputShift','outputShift');
% save('NNDataset_Raw_EPR.mat','dataInputApp','dataOutputApp','-append');
% save('NNDataset_P_EPR.mat','dataInputAppP','dataOutputAppP','-append');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Testing the trained ANN against each data set
% figure
% for i=1:N
%     disRate = diff(dataStruct(i).dis)./diff(dataStruct(i).time);
%     disRate = [disRate ; disRate(end)];
%     dis = dataStruct(i).dis;
%     neuralInput = [dataStruct(i).dis disRate];
%     subplot(2,4,i); plot(dataStruct(i).dis,dataStruct(i).load,...
%                         dataStruct(i).dis,NeuralNetwork_SR(neuralInput));
% %                         dataStruct(i).dis,SRNeuralNetwork_Single(strain));
%     clear disRate;
%     clear dis
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%