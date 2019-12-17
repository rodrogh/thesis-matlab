%%
%Script for training nueral networks of different topologies
clear all
clc
close all

%%
%The N datasets contain within a test must be concatenated vertically into
%a single dataset for training the NNs
load('TrainingSetProcStressIndV4.mat');

%% NNV4.0 
%The amount of data extracted from each test is fixed for all of and 
%dependent on the length of the smallest test (Raw Data)
%% Resize Rubber Bands Data (From 7 to 18)
%Resize is performed in compilePost.m directly
%% Variables Initialization
%Defining placeholder
trainData = PmatData;
validationData = PmatData;
%Defining number of materials
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR';'Nat100R'}
fields.N = fieldnames(PmatData);
%Find number of tests performance per type of test
for i=1:N
    %Defining number of type of tensile strength tests (3)
    M = length( fieldnames(PmatData. ( fields.N{i} ) ) );
    fields.M = fieldnames(PmatData.( fields.N{i} ) );
    for j=1:M
        testNum(j) = length( PmatData. ( fields.N{i} ).( fields.M{j}) );
        [~, tempIndex] = max(testNum);
        largestTest{i} = fields.M{tempIndex}; 
    end
end

for i=1:N
    M = length( fieldnames(PmatData. ( fields.N{i} ) ) );
    fields.M = fieldnames(PmatData.( fields.N{i} ) );
    trainData.(fields.N{i}) = [];
    validationData.(fields.N{i}) = [];
    concat_all = struct( [] );
    temp = struct( [] );    
    for j=1:M
        P = length( PmatData. ( fields.N{i} ).( fields.M{j}) );
        smallest = PmatData.(fields.N{i}).(fields.M{j})(1).stats.smallest;
        switch fields.M{j}
            case "disR50" 
                disR = 50/60; %milimeters per minute
            case "disR250" 
                disR = 250/60;
            case "disR500" 
                disR = 500/60;
        end
        strainR = disR/33; %Specimen Initial Length (mm / mm)
        for k=1:P
            %New Feature (Sep/19) extract same amount of data from all tests
            %Structure has different variables now:
            %pload is the processed data (I am using this to train the NN
            %first)
            %rload is the raw data from experiments
            %sload is the processed and smoothed data
            dataSize = length( PmatData.(fields.N{i}).(fields.M{j})(k).pload ); 
            rangeUnified = round( linspace(1 , dataSize, smallest ) );
            padding = 10;
            thickness = PmatData.(fields.N{i}).(fields.M{j})(k).stats.thickness;
            tempDisR = ones(1,smallest)*disR;
            tempStrainR = ones(1,smallest)*strainR;
            %% Load Dynamic Concatenation with padded zeros
            
            temp(j).thickness(k,:) = {thickness*ones(1,smallest+padding)};
            temp(j).load(k,:) = {padarray( PmatData.(fields.N{i}).(fields.M{j})(k).sload(rangeUnified), [0 padding],'pre' )};
            temp(j).dis(k,:) = {padarray( PmatData.(fields.N{i}).(fields.M{j})(k).pdis(rangeUnified), [0 padding],'pre')};
            temp(j).disR(k,:) = {padarray( tempDisR , [0 padding],'pre')};
            
            %% Stress Dynamic Concatenation with padded zeros
            temp(j).stress(k,:) = {padarray( PmatData.(fields.N{i}).(fields.M{j})(k).sstress(rangeUnified), [0 padding],'pre' )};
            temp(j).strain(k,:) = {padarray( PmatData.(fields.N{i}).(fields.M{j})(k).pstrain(rangeUnified), [0 padding],'pre' )};
            temp(j).strainR(k,:) = {padarray( tempStrainR , [0 padding],'pre')};
        end     
        %Before Horizontal concatenation: Extract the data from one test to
        %use it to measure the performance of the NN as an unknown dataset.
        %Since, the number of experiments differs from one strain rate to
        %the other, choose the strain rate with more tests to extract one
        %from there
        %The variable largestTest holds the experiment with the highest
        %number of tests
        if strcmp(fields.M{j},largestTest{i})
            %Horizontal concatenation
            %Validation Set (Unknown data)
            concat_all(j).uload = horzcat( temp(j).load{1} )';
            concat_all(j).udis = horzcat( temp(j).dis{1} )';
            concat_all(j).udisR = horzcat( temp(j).disR{1} )';
            concat_all(j).uthickness = horzcat( temp(j).thickness{1} )';
            concat_all(j).ustress = horzcat( temp(j).stress{1} )';
            concat_all(j).ustrain = horzcat( temp(j).strain{1} )';
            concat_all(j).ustrainR = horzcat( temp(j).strainR{1} )';
            %Training Set
            concat_all(j).load = horzcat( temp(j).load{2:end} )';
            concat_all(j).dis = horzcat( temp(j).dis{2:end} )';
            concat_all(j).disR = horzcat( temp(j).disR{2:end} )';
            concat_all(j).thickness = horzcat( temp(j).thickness{2:end} )';
            concat_all(j).stress = horzcat( temp(j).stress{2:end} )';
            concat_all(j).strain = horzcat( temp(j).strain{2:end} )';
            concat_all(j).strainR = horzcat( temp(j).strainR{2:end} )';
        else
            %Horizontal concatenation
            concat_all(j).load = horzcat( temp(j).load{:} )';
            concat_all(j).dis = horzcat( temp(j).dis{:} )';
            concat_all(j).disR = horzcat( temp(j).disR{:} )';
            concat_all(j).thickness = horzcat( temp(j).thickness{:} )';
            concat_all(j).stress = horzcat( temp(j).stress{:} )';
            concat_all(j).strain = horzcat( temp(j).strain{:} )';
            concat_all(j).strainR = horzcat( temp(j).strainR{:} )';
        end
        
    end
    
    %Extra concatenation of different displacement rates
    load_all = vertcat(concat_all(:).load);
    dis_all = vertcat(concat_all(:).dis);
    disR_all = vertcat(concat_all(:).disR);
    thickness_all = vertcat(concat_all(:).thickness);
    stress_all = vertcat(concat_all(:).stress);
    strain_all = vertcat(concat_all(:).strain);
    strainR_all = vertcat(concat_all(:).strainR);
    %
    uload_all = vertcat(concat_all(:).uload);
    udis_all = vertcat(concat_all(:).udis);
    udisR_all = vertcat(concat_all(:).udisR);
    uthickness_all = vertcat(concat_all(:).uthickness);
    ustress_all = vertcat(concat_all(:).ustress);
    ustrain_all = vertcat(concat_all(:).ustrain);
    ustrainR_all = vertcat(concat_all(:).ustrainR);
    
    %Shifted data (past/future values)
    load_shift1 = circshift(load_all,1,1);
    load_shift1(1) = 0;
    load_shift2 = circshift(load_all,2,1);
    load_shift2(1:2) = 0;
    
    dis_shift1 = circshift(dis_all,1,1);
    dis_shift1(1) = 0;
    dis_shift2 = circshift(dis_all,2,1);
    dis_shift2(1:2) = 0;
    
    disR_shift1 = circshift(disR_all,1,1);
    disR_shift1(1) = 0;
    disR_shift2 = circshift(disR_all,2,1);
    disR_shift2(1:2) = 0;
    
    stress_shift1 = circshift(stress_all,1,1);
    stress_shift1(1) = 0;
    stress_shift2 = circshift(stress_all,2,1);
    stress_shift2(1:2) = 0;
    
    strain_shift1 = circshift(strain_all,1,1);
    strain_shift1(1) = 0;
    strain_shift2 = circshift(strain_all,2,1);
    strain_shift2(1:2) = 0;
    
    strainR_shift1 = circshift(strainR_all,1,1);
    strainR_shift1(1) = 0;
    strainR_shift2 = circshift(strainR_all,2,1);
    strainR_shift2(1:2) = 0;
    
    
    %Shifted data (past/future values)
    uload_shift1 = circshift(uload_all,1,1);
    uload_shift1(1) = 0;
    uload_shift2 = circshift(uload_all,2,1);
    uload_shift2(1:2) = 0;
    
    udis_shift1 = circshift(udis_all,1,1);
    udis_shift1(1) = 0;
    udis_shift2 = circshift(udis_all,2,1);
    udis_shift2(1:2) = 0;
    
    udisR_shift1 = circshift(udisR_all,1,1);
    udisR_shift1(1) = 0;
    udisR_shift2 = circshift(udisR_all,2,1);
    udisR_shift2(1:2) = 0;
    
    ustress_shift1 = circshift(ustress_all,1,1);
    ustress_shift1(1) = 0;
    ustress_shift2 = circshift(ustress_all,2,1);
    ustress_shift2(1:2) = 0;
    
    ustrain_shift1 = circshift(ustrain_all,1,1);
    ustrain_shift1(1) = 0;
    ustrain_shift2 = circshift(ustrain_all,2,1);
    ustrain_shift2(1:2) = 0;
    
    ustrainR_shift1 = circshift(ustrainR_all,1,1);
    ustrainR_shift1(1) = 0;
    ustrainR_shift2 = circshift(ustrainR_all,2,1);
    ustrainR_shift2(1:2) = 0;
    
    %Saving data into structure
    trainData.(fields.N{i}).input_load = [dis_all dis_shift1 dis_shift2 ...
        disR_all disR_shift1 disR_shift2 thickness_all];
    trainData.(fields.N{i}).output_load = [load_all load_shift1 load_shift2];
    
    validationData.(fields.N{i}).input_load = [udis_all udis_shift1 udis_shift2 ...
        udisR_all udisR_shift1 udisR_shift2 uthickness_all];
    validationData.(fields.N{i}).output_load = [uload_all uload_shift1 uload_shift2];
    
    trainData.(fields.N{i}).input_stress = [strain_all strain_shift1 strain_shift2 ...
        strainR_all strainR_shift1 strainR_shift2 thickness_all];
    trainData.(fields.N{i}).output_stress = [stress_all stress_shift1 stress_shift2];
    
    validationData.(fields.N{i}).input_stress = [ustrain_all ustrain_shift1 ustrain_shift2 ...
        ustrainR_all ustrainR_shift1 ustrainR_shift2 uthickness_all];
    validationData.(fields.N{i}).output_stress = [ustress_all ustress_shift1 ustress_shift2];
end
%Save training data into file
save('FormattedTrainingSetIndV4.mat','trainData','validationData');
disp('Done')
