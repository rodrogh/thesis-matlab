%Create a .mat file containing a 1x3 variable with the postprocessed data
%from the .csv files in the form of: time, dis, load
clear all;
clc;
close all;
%Select which test data will be compiled
% test_to_compile = 'Stress Relaxation';
test_to_compile = 'Tensile Strength';
dataVersion = "V4";
matData = compileDataPaths(test_to_compile);

%% Manipulate Unified Data
PmatData = compileLoad(matData,test_to_compile);
PmatData = resizeRubberBandsData(PmatData);
switch test_to_compile
    case 'Tensile Strength'
        fileName = strcat("../Experimental Work/TrainingSetProc",dataVersion,".mat" );
        save(fileName,'PmatData');
    case 'Stress Relaxation'
        fileName = strcat("../Experimental Work/SRelSetProc",dataVersion,".mat" );
        save(fileName,'PmatData');
end
disp(strcat(test_to_compile," Data Saved"));
[~] = loadToStress(PmatData,test_to_compile,dataVersion);

%% Manipulate Formatted Data for NN Training
PmatData = compileLoadforNN(matData,test_to_compile);
PmatData = resizeRubberBandsData(PmatData);
switch test_to_compile
    case 'Tensile Strength'
        fileName = strcat("../Experimental Work/TrainingSetProcInd",dataVersion,".mat" );
        save(fileName,'PmatData');
    case 'Stress Relaxation'
        fileName = strcat("../Experimental Work/SRelSetProcInd",dataVersion,".mat" );
        save(fileName,'PmatData');
end
disp(strcat(test_to_compile," Data Saved"));
[~] = loadToStressforNN(PmatData,test_to_compile,dataVersion);
disp('Done');

function PmatData = compileLoad(matData,test_to_compile)
PmatData = matData;
N = length( fieldnames(matData) );
fields.N = fieldnames(matData);
%Function definition
%[pload, rload, rdis, rtime] = readInstronTable(N,iRow,iCol,paths,test)
for i=1:N
    disp(fields.N{i});
    M = length( fieldnames(matData. ( fields.N{i} ) ) );
    fields.M = fieldnames(matData.(fields.N{i}));
    for j=1:M
        disp(fields.M{j})
        P = length( matData.(fields.N{i}).(fields.M{j}).paths );
        %Slightly modify matData layout to include pload, rload, rdis, rtime
        %outputs of the readInstronTable
        [ PmatData.(fields.N{i}).(fields.M{j}).rload ,...
            PmatData.(fields.N{i}).(fields.M{j}).rdis ,...
            PmatData.(fields.N{i}).(fields.M{j}).rtime ,...
            PmatData.(fields.N{i}).(fields.M{j}).pload ,...
            PmatData.(fields.N{i}).(fields.M{j}).pdis ,...
            PmatData.(fields.N{i}).(fields.M{j}).ptime,...
            PmatData.(fields.N{i}).(fields.M{j}).sload,... %smoothed load
            PmatData.(fields.N{i}).(fields.M{j}).stats ] = ...%Stats Values
            readInstronTable(P,2,0 ,...
            matData.(fields.N{i}).(fields.M{j}).paths ,...
            test_to_compile);
        if i==5 %Thickness of PR is 6
            PmatData.(fields.N{i}).(fields.M{j}).stats.thickness = 6;  %Defaul Thickness
        else
            PmatData.(fields.N{i}).(fields.M{j}).stats.thickness = 1.5;  %Defaul Thickness
        end
    end
end

end

function PmatData = compileLoadforNN(matData,test_to_compile)
PmatData = matData;
N = length( fieldnames(matData) );
fields.N = fieldnames(matData);
%Function definition
%[pload, rload, rdis, rtime] = readInstronTable(N,iRow,iCol,paths,test)
for i=1:N
    M = length( fieldnames(matData. ( fields.N{i} ) ) );
    fields.M = fieldnames(matData.(fields.N{i}));
    for j=1:M
        P = length( matData.(fields.N{i}).(fields.M{j}).paths );
        %The k loop is only for Neural Networks. The aim is to smooht all
        %the available raw data individually and use it for training. When
        %this is not required, pass the variable path to the
        %readInstronTable(N,iRow,iCol,paths,test) function as usual create
        %a smoothed dataset per experimental setup.
        paths = matData.(fields.N{i}).(fields.M{j}).paths;
        u = zeros(1,P);
        for k=1:P
            %Overwrite the Nx1 cells path variable to comply with the
            %other parameters
            PmatData.(fields.N{i}).(fields.M{j})(k).paths = paths{k};
            [ PmatData.(fields.N{i}).(fields.M{j})(k).rload ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).rdis ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).rtime ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).pload ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).pdis ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).ptime ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).sload ,...
                PmatData.(fields.N{i}).(fields.M{j})(k).stats] = ... %Stats values
                readInstronTable(1,2,0 ,...
                paths(k) ,...
                test_to_compile);
            %Add relevant values to stats structure
            if i==5 %Thickness of PR is 6
                PmatData.(fields.N{i}).(fields.M{j})(k).stats.thickness = 6;  %Defaul Thickness
            else
                PmatData.(fields.N{i}).(fields.M{j})(k).stats.thickness = 1.5;  %Defaul Thickness
            end
            u(k) = size(PmatData.(fields.N{i}).(fields.M{j})(k).pload,2);
        end
        smallest = min(u);
        for k=1:P
            PmatData.(fields.N{i}).(fields.M{j})(k).stats.smallest = smallest;
        end
        
    end
end
end

function PmatData = loadToStress(matData,test_to_compile,dataVersion)
PmatData = matData;
specimenLength = 33; %mm
specimenWidth = 6; %mm
N = length( fieldnames(matData) );
fields.N = fieldnames(matData);
%Function definition
%[pload, rload, rdis, rtime] = readInstronTable(N,iRow,iCol,paths,test)
for i=1:N
    M = length( fieldnames(matData. ( fields.N{i} ) ) );
    fields.M = fieldnames(matData.(fields.N{i}));
    for j=1:M
        P = length( matData.(fields.N{i}).(fields.M{j}) );
        %% IMPORTANT: LOAD AND DISPLACEMENT DATA FOR RUBBER BANDS IN THIS
        %% FORMAT ARE NOT RELEVANT. ONLY THE STRESS DATA IS BECAUSE IT TAKES
        %% INTO ACCOUNT THE MATERIALS DIMENSIONS
        %Test from different rubber bands have different lengths,
        %find them
        for k=1:P
            u1(k) = length(matData.(fields.N{i}).(fields.M{j})(k).rload);
            u2(k) = length(matData.(fields.N{i}).(fields.M{j})(k).pload);
            u3(k) = length(matData.(fields.N{i}).(fields.M{j})(k).sload);
        end
        %Sample the same data from each test, based on the smallest
        %test
        for k=1:P
            range1 = round(linspace(1,u1(k),min(u1)));
            range2 = round(linspace(1,u2(k),min(u2)));
            range3 = round(linspace(1,u3(k),min(u3)));
            % Placeholder to improve readability of code
            thickness = matData.(fields.N{i}).(fields.M{j})(k).stats.thickness;
            crossArea = (specimenWidth*thickness*1e-6);
            % Compile data and calculate stress
            temp_rstress(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).rload(range1)...
                /crossArea;
            temp_rstrain(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).rdis(range1)...
                / specimenLength;
            temp_rtime(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).rtime(range1);
            temp_pstress(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).pload(range2)...
                /crossArea;
            temp_pstrain(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).pdis(range2)...
                / specimenLength;
            temp_ptime(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).ptime(range2);
            temp_sstress(k,:) = ...
                matData.(fields.N{i}).(fields.M{j})(k).sload(range3)...
                /crossArea;
            % Extra parameters
            temp_eng_stress = matData.(fields.N{i}).(fields.M{j})(k).stats.medLoad...
                /crossArea;
            temp_eng_strain = matData.(fields.N{i}).(fields.M{j})(k).stats.medDis...
                /specimenLength;
            temp_eng_stress_end = matData.(fields.N{i}).(fields.M{j})(k).stats.minLoad...
                /crossArea;
        end
        PmatData.(fields.N{i}).(fields.M{j}) = [];
        %Store Averaged Stress values of all tests in all variables of
        %the structure to retain consistency (duplicating)
        PmatData.(fields.N{i}).(fields.M{j}).rstress = ...
            mean( temp_rstress,1  );
        PmatData.(fields.N{i}).(fields.M{j}).pstress = ...
            mean( temp_pstress,1  );
        PmatData.(fields.N{i}).(fields.M{j}).sstress  = ...
            mean( temp_sstress,1  );
        PmatData.(fields.N{i}).(fields.M{j}).rstrain = ...
            mean( temp_rstrain,1  );
        PmatData.(fields.N{i}).(fields.M{j}).pstrain = ...
            mean( temp_pstrain,1  );
        PmatData.(fields.N{i}).(fields.M{j}).rtime = ...
            mean( temp_rtime,1  );
        PmatData.(fields.N{i}).(fields.M{j}).ptime = ...
            mean( temp_ptime,1  );
        %Depending on the test, this value has two different meanings
        %Tensile Strength - Engineering Stress and Strain
        %Stress Relaxation - Initial Stress (median)
        PmatData.(fields.N{i}).(fields.M{j}).eng_stress = ...
            mean( temp_eng_stress, 1 );
        PmatData.(fields.N{i}).(fields.M{j}).eng_strain = ...
            mean( temp_eng_strain, 1 );        
        %Stress Relaxation - Stress at the end of the test (median)
        PmatData.(fields.N{i}).(fields.M{j}).eng_stress_end = ...
            mean( temp_eng_stress_end, 1 );
        
        clear temp_rstress temp_rstrain ...
            temp_rtime temp_ptime ...
            temp_pstress temp_pstrain temp_sstress ...
            u1 u2 u3
    end
end
switch test_to_compile
    case 'Tensile Strength'
        fileName = strcat("../Experimental Work/TrainingSetProcTrainingSetProcStress",dataVersion,".mat" );
        save(fileName,'PmatData');
    case 'Stress Relaxation'
        fileName = strcat("../Experimental Work/TrainingSetProcSRelSetProcStress",dataVersion,".mat" );
        save(fileName,'PmatData');
end
disp(strcat(test_to_compile," Data Saved"));
end

function PmatData = loadToStressforNN(matData,test_to_compile,dataVersion)
PmatData = matData;
specimenLength = 33; %mm
specimenWidth = 6; %mm
N = length( fieldnames(matData) );
fields.N = fieldnames(matData);
%Function definition
%[pload, rload, rdis, rtime] = readInstronTable(N,iRow,iCol,paths,test)
for i=1:N
    M = length( fieldnames(matData. ( fields.N{i} ) ) );
    fields.M = fieldnames(matData.(fields.N{i}));
    for j=1:M
        P = length( matData.(fields.N{i}).(fields.M{j}) );
        for k=1:P
            %Extract material thickness
            thickness = matData.(fields.N{i}).(fields.M{j})(k).stats.thickness;            
            PmatData.(fields.N{i}).(fields.M{j})(k).rstress = ...
                matData.(fields.N{i}).(fields.M{j})(k).rload / ...
                (specimenWidth*thickness*1e-6); %Pascals
            PmatData.(fields.N{i}).(fields.M{j})(k).pstress = ...
                matData.(fields.N{i}).(fields.M{j})(k).pload / ...
                (specimenWidth*thickness*1e-6); %Pascals
            PmatData.(fields.N{i}).(fields.M{j})(k).sstress = ...
                matData.(fields.N{i}).(fields.M{j})(k).sload / ...
                (specimenWidth*thickness*1e-6); %Pascals            
            PmatData.(fields.N{i}).(fields.M{j})(k).rstrain = ...
                matData.(fields.N{i}).(fields.M{j})(k).rdis / specimenLength; %Not Percentage
            PmatData.(fields.N{i}).(fields.M{j})(k).pstrain = ...
                matData.(fields.N{i}).(fields.M{j})(k).pdis / specimenLength; %Not Percentage
        end
    end
end
switch test_to_compile
    case 'Tensile Strength'
        fileName = strcat("../Experimental Work/TrainingSetProcTrainingSetProcStressInd",dataVersion,".mat" );
        save(fileName,'PmatData');
    case 'Stress Relaxation'
        fileName = strcat("../Experimental Work/TrainingSetProcSRelSetProcStressInd",dataVersion,".mat" );
        save(fileName,'PmatData');
end
disp(strcat(test_to_compile," Data Saved"));
end