function [matData] = compileDataPaths(test_name)
%Function to compile all the datasets of the available soft materials,
%grouping them per displacement rate (tensile strength).
%Output:
%   matData - structure containing the three diferent displacement rates
%   per tensile strength experiment (50,250,500) for all the six sof
%   materials (EPR,FR,NatR,NR,PR,SR)
%
%Optionally this function also saves the structure into DatasetPathsAll.mat

switch test_name
    case 'Stress Relaxation'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SR
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

SR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EPR

temp = {'ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO EPR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

EPR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FR

temp = {'ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO FR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

FR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NatR

temp = {'ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO NatR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

NatR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NR

temp = {'ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO NR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

NR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PR

temp = {'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

PR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBYellowThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBYellowThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBYellowThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

RBYellowThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBRedThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBRedThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBRedThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

RBRedThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBBlueThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBBlueThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBBlueThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

RBBlueThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBGreenThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBGreenThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBGreenThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

RBGreenThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBBlackThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBBlackThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBBlackThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

RBBlackThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBOrangeThin
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'}; %Thin
L180min.paths(1:length(temp),:) = temp;

RBOrangeThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RBOrangeThick
temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
L5min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
L15min.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'}; %Thick
L180min.paths(1:length(temp),:) = temp;

RBOrangeThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Tensile Strength'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SR
temp = {'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_11';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_12';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_13';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_14';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_15'};
disR50.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;

SR = struct('disR50' , {disR50}, 'disR250', {disR250});
% SR = {disR50 , disR250 , disR500};
clear disR50 disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EPR

temp = {'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_11';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_12';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_13';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_14';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_15';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_16'};
disR50.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_5'};
disR500.paths(1:length(temp),:) = temp;

EPR = struct('disR50' , {disR50} , 'disR500', {disR500});
% EPR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FR

temp = {'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_8'};
disR50.paths(1:length(temp),:) = temp;


temp = {  'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_8'};
disR250.paths(1:length(temp),:) = temp;


temp = {'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_5'};
disR500.paths(1:length(temp),:) = temp;

FR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% FR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NatR

temp = {'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_11'};
disR50.paths(1:length(temp),:) = temp;

temp = {   'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_5'};
disR250.paths(1:length(temp),:) = temp;


temp = {    'ASTM D412\Tensile\NatR\RDSO NatR Tensile 500.is_tens_RawData\Specimen_RawData_1'};
disR500.paths(1:length(temp),:) = temp;

NatR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% NatR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NR

temp = {'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_8'};
disR50.paths(1:length(temp),:) = temp;


temp = {    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;


temp = {'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_6'};
disR500.paths(1:length(temp),:) = temp;

NR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% NR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PR
temp = {'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_8';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_9';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_10';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_11';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_12';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_13'};
disR50.paths(1:length(temp),:) = temp;

temp = { 'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;


temp = { 'ASTM D412\Tensile\PR\RDSO PR6 Tensile 500.is_tens_RawData\Specimen_RawData_1' };
disR500.paths(1:length(temp),:) = temp;

PR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% PR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBYellowThin
temp = {'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 50.is_tens_RawData\Specimen_RawData_1'};
disR50.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_6'};
disR250.paths(1:length(temp),:) = temp;

RBYellowThin = struct('disR50' , {disR50}, 'disR250', {disR250} );
clear disR50 disR250 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBYellowThick

temp = {    'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_5'};
disR250.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 500.is_tens_RawData\Specimen_RawData_1'};
disR500.paths(1:length(temp),:) = temp;

RBYellowThick = struct( 'disR250', {disR250} , 'disR500', {disR500});
clear disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBRedThin
temp = {'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_2';   
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_8'};
disR250.paths(1:length(temp),:) = temp;

RBRedThin = struct( 'disR250', {disR250} );
clear disR250 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBRedThick
temp = {'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_3';    
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;

RBRedThick = struct( 'disR250', {disR250} );
clear disR250 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlueThin

temp = {'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_8'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlueThin = struct( 'disR250', {disR250} );
clear disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlueThick

temp = {'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_3';    
    'ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlueThick = struct( 'disR250', {disR250} );
clear disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBGreenThin

temp = {'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_5';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_6';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_7';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_8'};
    
disR500.paths(1:length(temp),:) = temp;

RBGreenThin = struct( 'disR500', {disR500} );
clear disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBGreenThick

temp = {'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_4'};
    
disR500.paths(1:length(temp),:) = temp;

RBGreenThick = struct( 'disR500', {disR500} );
clear disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlackThin
temp = {'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_6'};
    
disR250.paths(1:length(temp),:) = temp;

temp = {'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 500.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 500.is_tens_RawData\Specimen_RawData_2'};
    
disR500.paths(1:length(temp),:) = temp;

RBBlackThin = struct( 'disR250', {disR250} , 'disR500', {disR500});
clear disR250 disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlackThick
temp = {'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_1';
    'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlackThick = struct( 'disR250', {disR250} );
clear disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBOrangeThin

temp = {'ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_2';
    'ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_3';
    'ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_4';
    'ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBOrangeThin = struct( 'disR250', {disR250} );
clear disR250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBOrangeThick

temp = {'ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_1'};
    
disR250.paths(1:length(temp),:) = temp;

RBOrangeThick = struct( 'disR250', {disR250} );
clear disR250;
end

matData = struct('EPR',EPR,'FR',FR,'NatR',NatR,'NR',NR,'PR',PR,'SR',SR,...
    'RBYellowThin',RBYellowThin, 'RBYellowThick',RBYellowThick,...
    'RBRedThin',RBRedThin, 'RBRedThick',RBRedThick,...
    'RBBlueThin',RBBlueThin , 'RBBlueThick',RBBlueThick ,...
    'RBGreenThin',RBGreenThin,'RBGreenThick',RBGreenThick,...
    'RBBlackThin',RBBlackThin,'RBBlackThick',RBBlackThick,...
    'RBOrangeThin',RBOrangeThin,'RBOrangeThick',RBOrangeThick);

% save('DatasetPathsAll.mat','matData');
end

