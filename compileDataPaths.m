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
temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

SR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EPR

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO EPR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

EPR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FR

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO FR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

FR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NatR

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NatR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

NatR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NR

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO NR SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

NR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PR

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 50.is_trelax_RawData\Specimen_RawData_5'};
L5min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_2';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_3';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_4';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 100.is_trelax_RawData\Specimen_RawData_5'};
L15min.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\PR\RDSO PR6 SRelax 500.is_trelax_RawData\Specimen_RawData_1'};
L180min.paths(1:length(temp),:) = temp;

PR = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
clear L5min L15min L180min;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resistance Bands - Compile the available data into the same set from here
temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBBlack = struct('L15min',{L180min});
clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1';
    '../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
L180min.paths(1:length(temp),:) = temp;

RBOrange = struct('L15min',{L180min});
clear L5min L15min L180min;

matData = struct('EPR',EPR,'FR',FR,'NatR',NatR,'NR',NR,'PR',PR,'SR',SR,...
    'RBBlack',RBBlack,'RBOrange',RBOrange);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBYellowThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBYellowThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBYellowThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBYellowThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBRedThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBRedThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBRedThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBRedThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBBlueThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBBlueThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBBlueThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBBlueThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBGreenThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBGreenThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBGreenThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBGreenThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBBlackThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBBlackThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBBlackThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBBlack SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'};
% L180min.paths(1:length(temp),:) = temp;
% 
% RBBlackThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBOrangeThin
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_2'}; %Thin
% L180min.paths(1:length(temp),:) = temp;
% 
% RBOrangeThin = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % RBOrangeThick
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 50.is_trelax_RawData\Specimen_RawData_1'};
% L5min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO SR SRelax 100.is_trelax_RawData\Specimen_RawData_1'};
% L15min.paths(1:length(temp),:) = temp;
% 
% temp = {'../Experimental Work/../Experimental Work/ASTM D412\Stress Relaxation\RDSO RBOrange SRel 10d 1000s.is_trelax_RawData\Specimen_RawData_1'}; %Thick
% L180min.paths(1:length(temp),:) = temp;
% 
% RBOrangeThick = struct('L5min',{L5min},'L15min',{L15min},'L180min',{L180min});
% clear L5min L15min L180min;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Tensile Strength'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SR
temp = {'../Experimental Work/../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_11';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_12';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_13';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_14';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 50.is_tens_RawData\Specimen_RawData_15'};
disR50.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\SR\RDSO SR Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;

SR = struct('disR50' , {disR50}, 'disR250', {disR250});
% SR = {disR50 , disR250 , disR500};
clear disR50 disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EPR

temp = {'../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_11';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_12';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_13';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_14';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_15';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 50.is_tens_RawData\Specimen_RawData_16'};
disR50.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\EPR\RDSO EPR Tensile 500.is_tens_RawData\Specimen_RawData_5'};
disR500.paths(1:length(temp),:) = temp;

EPR = struct('disR50' , {disR50} , 'disR500', {disR500});
% EPR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%FR

temp = {'../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 50.is_tens_RawData\Specimen_RawData_8'};
disR50.paths(1:length(temp),:) = temp;


temp = {  '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 250.is_tens_RawData\Specimen_RawData_8'};
disR250.paths(1:length(temp),:) = temp;


temp = {'../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\FR\RDSO FR Tensile 500.is_tens_RawData\Specimen_RawData_5'};
disR500.paths(1:length(temp),:) = temp;

FR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% FR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NatR

temp = {'../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_8';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_9';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_10';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 50.is_tens_RawData\Specimen_RawData_11'};
disR50.paths(1:length(temp),:) = temp;

temp = {   '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 250.is_tens_RawData\Specimen_RawData_5'};
disR250.paths(1:length(temp),:) = temp;


temp = {    '../Experimental Work/ASTM D412\Tensile\NatR\RDSO NatR Tensile 500.is_tens_RawData\Specimen_RawData_1'};
disR500.paths(1:length(temp),:) = temp;

NatR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% NatR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%NR

temp = {'../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 50.is_tens_RawData\Specimen_RawData_8'};
disR50.paths(1:length(temp),:) = temp;


temp = {    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;


temp = {'../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\NR\RDSO NR Tensile 500.is_tens_RawData\Specimen_RawData_6'};
disR500.paths(1:length(temp),:) = temp;

NR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% NR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PR
temp = {'../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_8';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_9';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_10';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_11';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_12';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 50.is_tens_RawData\Specimen_RawData_13'};
disR50.paths(1:length(temp),:) = temp;

temp = { '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;


temp = { '../Experimental Work/ASTM D412\Tensile\PR\RDSO PR6 Tensile 500.is_tens_RawData\Specimen_RawData_1' };
disR500.paths(1:length(temp),:) = temp;

PR = struct('disR50' , {disR50}, 'disR250', {disR250} , 'disR500', {disR500});
% PR = {disR50 , disR250 , disR500};
clear disR50 disR250 disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBYellowThin
temp = {'../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 50.is_tens_RawData\Specimen_RawData_1'};
disR50.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_6'};
disR250.paths(1:length(temp),:) = temp;

RBYellowThin = struct('disR50' , {disR50}, 'disR250', {disR250} );
clear disR50 disR250 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBYellowThick

temp = {    '../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 250.is_tens_RawData\Specimen_RawData_5'};
disR250.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/ASTM D412\Tensile\RBYellow\RDSO RBYellow Tensile 500.is_tens_RawData\Specimen_RawData_1'};
disR500.paths(1:length(temp),:) = temp;

RBYellowThick = struct( 'disR250', {disR250} , 'disR500', {disR500});
clear disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBRedThin
temp = {'../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_2';   
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_8'};
disR250.paths(1:length(temp),:) = temp;

RBRedThin = struct( 'disR250', {disR250} );
clear disR250 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBRedThick
temp = {'../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_3';    
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\RBRed\RDSO RBRed Tensile 250.is_tens_RawData\Specimen_RawData_7'};
disR250.paths(1:length(temp),:) = temp;

RBRedThick = struct( 'disR250', {disR250} );
clear disR250 ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlueThin

temp = {'../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_8'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlueThin = struct( 'disR250', {disR250} );
clear disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlueThick

temp = {'../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_3';    
    '../Experimental Work/ASTM D412\Tensile\RBBlue\RDSO RBBlue Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlueThick = struct( 'disR250', {disR250} );
clear disR250;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBGreenThin

temp = {'../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_5';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_6';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_7';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_8'};
    
disR500.paths(1:length(temp),:) = temp;

RBGreenThin = struct( 'disR500', {disR500} );
clear disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBGreenThick

temp = {'../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\RBGreen\RDSO RBGreen Tensile 500.is_tens_RawData\Specimen_RawData_4'};
    
disR500.paths(1:length(temp),:) = temp;

RBGreenThick = struct( 'disR500', {disR500} );
clear disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlackThin
temp = {'../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_6'};
    
disR250.paths(1:length(temp),:) = temp;

temp = {'../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 500.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 500.is_tens_RawData\Specimen_RawData_2'};
    
disR500.paths(1:length(temp),:) = temp;

RBBlackThin = struct( 'disR250', {disR250} , 'disR500', {disR500});
clear disR250 disR500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBBlackThick
temp = {'../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_1';
    '../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\RBBlack\RDSO RBBlack Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBBlackThick = struct( 'disR250', {disR250} );
clear disR250 disR500;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBOrangeThin

temp = {'../Experimental Work/ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_2';
    '../Experimental Work/ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_3';
    '../Experimental Work/ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_4';
    '../Experimental Work/ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_5'};
    
disR250.paths(1:length(temp),:) = temp;

RBOrangeThin = struct( 'disR250', {disR250} );
clear disR250;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%RBOrangeThick

temp = {'../Experimental Work/ASTM D412\Tensile\RBOrange\RDSO RBOrange Tensile 250.is_tens_RawData\Specimen_RawData_1'};
    
disR250.paths(1:length(temp),:) = temp;

RBOrangeThick = struct( 'disR250', {disR250} );
clear disR250;

matData = struct('EPR',EPR,'FR',FR,'NatR',NatR,'NR',NR,'PR',PR,'SR',SR,...
    'RBYellowThin',RBYellowThin, 'RBYellowThick',RBYellowThick,...
    'RBRedThin',RBRedThin, 'RBRedThick',RBRedThick,...
    'RBBlueThin',RBBlueThin , 'RBBlueThick',RBBlueThick ,...
    'RBGreenThin',RBGreenThin,'RBGreenThick',RBGreenThick,...
    'RBBlackThin',RBBlackThin,'RBBlackThick',RBBlackThick,...
    'RBOrangeThin',RBOrangeThin,'RBOrangeThick',RBOrangeThick);

end

% save('DatasetPathsAll.mat','matData');
end

