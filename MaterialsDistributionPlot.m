clear all;
close all
clc;
%Script to plot a curve of Ultimate Load/displacement of all the materials
load('TrainingSetProc.mat');

%Defining number of materials
N = length( fieldnames(PmatData) );
%Materials NAMES: {'EPR';'FR';'NatR';'NR';'PR';'SR'}
fields.N = fieldnames(PmatData);
%Statistic values extracted from:
fields.M = {'disR500';'disR500';'disR250';...
            'disR500';'disR250';'disR250';...
            'disR250';'disR250';...
            'disR250';'disR250';...
            'disR250';'disR250';...
            'disR500';'disR500';...
            'disR250';'disR250';...
            'disR250';'disR250'};
for i=1:N
    x = PmatData.(fields.N{i}).(fields.M{i}).stats.minDis;
    y = PmatData.(fields.N{i}).(fields.M{i}).stats.minLoad;
    w = PmatData.(fields.N{i}).(fields.M{i}).stats.maxDis - x;
    h = PmatData.(fields.N{i}).(fields.M{i}).stats.maxLoad - y;
    color = rand(1,3);
    h = rectangle('Position', [x y w h], 'FaceColor', color ,'Curvature', [1 1]);
    hline(i) = line(NaN,NaN,'LineWidth',2,'Color',color);
end
legend(hline,fields.N);