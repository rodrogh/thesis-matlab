%Radial Basis Function Test
%Parameters as used in the literature
%Data is not normalized

clear all;
close all;
clc;
%% Initialization: Loading data and preparing assets
load('FormattedTrainingSetIndV3.mat'); %trainData
fields.N = fieldnames(trainData); %Field names (materials tested)
N = length( fields.N ); %Number of materials tested
%% Network Design: Data division in 80/10/10
%Focus on Natural Rubber
temp = trainData.NatR.inputs(:,1)';
Q = size(temp, 2); %Number of samples
Q1 = floor(Q * 0.80);
Q2 = floor( (Q - Q1)/2 );
Q3 = Q - Q1 - Q2;
ind = randperm(Q);
ind1 = ind(1:Q1);
ind2 = ind(Q1 + (1:Q2));
ind3 = ind(Q1 + Q2 + (1:Q3));
%% Network Design: Defining Topologies to Test(I/O)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Inputs: displacement and displacement rate
%Outputs: load
input_set = vertcat(trainData.NatR.inputs(:,1)',trainData.NatR.inputs(:,4)');
output_set = vertcat(trainData.NatR.outputs(:,1)');
input_div = {input_set(:,ind1),... %Training set
                input_set(:,ind2),... %Validation set
                input_set(:,ind3) };   %Test set
output_div = {output_set(:,ind1),... %Training set
                output_set(:,ind2),... %Validation set
                output_set(:,ind3) };   %Test set
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Network Design: Parameters for Training of Radial Basis Function Network
goal = 0.1; %literature
spread = 7; %literature
mn = 100; %Maximum number of neurons (default)
df = 10;%Number of neurons to add between displays (default)
%% Network Training
net = newrb( [input_div{1,1} input_div{1,2}],...
            [output_div{1,1} output_div{1,2}],...
            goal,spread,mn,df);
%% Network Testing
out_p = sim( net , input_div{1,3} );
perf = perform(net,output_div{1,3},out_p);

function BayesianRegTrain(trainData,matDesired,matNames)
%% Initialization

%% Network Design: Inputs and Outputs
% inputs = 
% outputs =
%% Network Design

end

