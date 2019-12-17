%% Description:
%Scripts to identify the toe and elastic regions of all the soft materials.
clear all
clc
close all
%Strings for plotting
f_names.mats_name = {'Ethylene Polypropylene Rubber';
    'Fluorocarbon Rubber';
    'Natural Rubber with Polyester';
    'Nitrile Rubber';
    'Polyethylene Rubber';
    'Silicone Rubber';
    'Natural Rubber 100%'};

%% Calculating ultimate values and yield offset
[tensileStrength] = getTensileStrength();
[tensileStrengthInd] = getTensileStrengthBandsInd();
%% Ploting Stress-Strain Curves
% plotTensileStrength(tensileStrength);
%% Ploting Offset Yield Strength
% plotOffsetYieldStrength(tensileStrength,f_names);
%% Export TensileStrength Properties to a table
% tensileTable = getFormattedTable(tensileStrength);
% tensileTableInd = getFormattedTableInd(tensileStrengthInd);
[X,Y,Z]=getThicknessSurface(tensileStrengthInd);
disp('Done');

function [X,Y,Z] = getThicknessSurface(tensileStrength)
% There only viable solution is to look at scatterInterpolat. Mapping the
% existing values in a larger space will cause duplicates and hence abrut
% changes in the plot. Interpolation could deal with this. Put a pin on
% this task, it is not critical.
testLengths = [];
maxStress = [];
minStress = [];
maxStrain = [];
minStrain = [];
thickness = [];
%% Extract the dimensions of each test and then sample the same amount of dat from each test
fields.N = fieldnames(tensileStrength);
N = length(fields.N);
for i=N %Only for rubber bands
    fields.M = fieldnames(tensileStrength.(fields.N{i}));
    M = length(fields.M);
    for j=2
        fields.P = fieldnames(tensileStrength.(fields.N{i}).(fields.M{j}));
        P = length(tensileStrength.(fields.N{i}).(fields.M{j}));
        for k=1:P
            testLengths(end+1) = ...
                length( tensileStrength.(fields.N{i}).(fields.M{j})(k).sstress);
            maxStress(end+1) = max(tensileStrength.(fields.N{i}).(fields.M{j})(k).sstress);
            minStress(end+1) = min(tensileStrength.(fields.N{i}).(fields.M{j})(k).sstress);
            maxStrain(end+1) = max(tensileStrength.(fields.N{i}).(fields.M{j})(k).pstrain);
            minStrain(end+1) = min(tensileStrength.(fields.N{i}).(fields.M{j})(k).pstrain);
            thickness(end+1) = tensileStrength.(fields.N{i}).(fields.M{j})(k).thickness;
        end
    end
end
%% Sample same amount of data from each test
index = 1;
% sampleSize = min(testLengths);
sampleSize = 25;
meshSize = 10*sampleSize;
[X,Y] = meshgrid( linspace( min(minStrain), max(maxStrain),meshSize ) ,...
    linspace( min(thickness), max(thickness),meshSize ));
strain = zeros(length(testLengths),sampleSize );
stress = zeros(length(testLengths),sampleSize );
for i=N %Only for rubber bands
    fields.M = fieldnames(tensileStrength.(fields.N{i}));
    M = length(fields.M);
    for j=2 %Consider only the test where you have all the thickness
        fields.P = fieldnames(tensileStrength.(fields.N{i}).(fields.M{j}));
        P = length(tensileStrength.(fields.N{i}).(fields.M{j}));
        for k=1:P
            sampleIndex = round( linspace(1, testLengths(index) , sampleSize ) );
            strain(index,:) = tensileStrength.(fields.N{i}).(fields.M{j})(k).pstrain(sampleIndex);
            stress(index,:) = tensileStrength.(fields.N{i}).(fields.M{j})(k).sstress(sampleIndex);
            
            index = index +1;
        end
    end
end
%Sorting rows in ascending order depending on thickness
[Y,I] = sort(thickness); %The variable I holds the sorted indexes
X = strain(I,:);
Z = stress(I,:);
%I have to sort the rows depending on the band thickness
mesh(Z,X);
end

function tbl = getFormattedTable(tensileStrength)
fields.N = fieldnames(tensileStrength);
N = length(fields.N);
T_data = table();
T_test = table();
for i=1:N
    fields.M = fieldnames(tensileStrength.(fields.N{i}));
    M = length(fields.M);
    for j=1:M
        T_test = [T_test; fields.N(i),fields.M(j)];
        tempT = struct2table(tensileStrength.(fields.N{i}).(fields.M{j}),'AsArray',true);
        T_data = [ T_data ; tempT(1,6:10) ];
    end
end
T_main = [T_test,T_data];
T_main.Properties.VariableNames(1)={'Materials'};
T_main.Properties.VariableNames(2)={'StrainRate'};
tbl = T_main;
end

function tbl = getFormattedTableInd(tensileStrength)
fields.N = fieldnames(tensileStrength);
N = length(fields.N);
T_data = table();
T_test = table();
for i=N
    fields.M = fieldnames(tensileStrength.(fields.N{i}));
    M = length(fields.M);
    for j=1:M
        fields.P = fieldnames(tensileStrength.(fields.N{i}).(fields.M{j}));
        P = length(tensileStrength.(fields.N{i}).(fields.M{j}));
        for k=1:P
            matName = tensileStrength.(fields.N{i}).(fields.M{j})(k).name;
            T_test = [T_test; matName , fields.M(j)];
            tempT = struct2table(tensileStrength.(fields.N{i}).(fields.M{j})(k),'AsArray',true);
            T_data = [ T_data ; tempT(1,6:10) ];
        end
    end
end
T_main = [T_test,T_data];
T_main.Properties.VariableNames(1)={'Materials'};
T_main.Properties.VariableNames(2)={'StrainRate'};
tbl = T_main;
end

function [yield_off, yfit_el, yfit_off, el_mod] = offsetYield(x,y,offset,size_el)
%Apply linear regression in the first X% strain
s_range = find( x <= size_el);
s_range_length = length(s_range);
%Approach using linear fit (origin)
% dlm00 = fitlm(x(s_range),y(s_range),'Intercept',false); %Linear fit model passing through origin
% m = dlm00.Coefficients{1,1};
% yfit_el = x(1:s_range_length*2) * m; %Linear fit passing through origin
% diff = y - (x(:) - (offset*lo/100) ) * m; %Linear fit with offset
%Approach using polyfit
p = polyfit(x(s_range),y(s_range),1);
m = p(1);
yfit_el = polyval(p,x(1:s_range_length*2));
diff = y - (x(:) - offset ) * m; %Linear fit with offset
[~ , index] = min( abs(diff) );
%Coordinates (x,y) of the offset yield strength
yield_off = [x(index), y(index)];
%Linear fit with offset up to the intercept
yfit_off = ( x( 1:index ) - offset ) * m;
el_mod = m; %Modulus of elasticity
end

function tensileStrength = getTensileStrength()
%% Initialization
%Loading the tensile strength processed data. The file: TrainingSetProc.mat
%contains the averaged data of all the specimens used n a single type of
%tensile strenght test.
load('TrainingSetProcStressV4.mat');%Stress in Pa, and Strain NOT in %
%Defining placeholders
f_names.mats =  fieldnames(PmatData); %contain the names of all the fields in PmatData
mats = length( f_names.mats ); %materials
tensileStrength = PmatData; %will store tensile data per test
for i = 1:mats
    %Different tests were performed to each material. In here I find out
    %how many tests and their name.
    f_names.tests =  fieldnames(PmatData.(f_names.mats{i}));
    tests = length(  f_names.tests );
    for j = 1:tests
        %Inside the struct of each test there are many variables, as
        %follows: paths,pload,rload,rdis,rtime,stats and stress values
        f_names.data = fieldnames( PmatData.( f_names.mats{i} ).( f_names.tests{j} ));
        %Clear previous data to store tensile properties in placeholder
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}) = [];
        %% Ultimate Values
        %The averaged values for rubber bands are stored in different
        %variables        
        if i==7 %Rubber Bands Index
            rstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).mrstress'; %rstress
            pstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).mpstress'; %pstress
            sstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).msstress'; %sstress
            rstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).mrstrain'; %rstrain
            pstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).mpstrain'; %pstrain    
        else
            rstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstress'; %rstress
            pstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstress'; %pstress
            sstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).sstress'; %sstress
            rstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstrain'; %rstrain
            pstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstrain'; %pstrain
        end
        %Calculate Ultimate Values
        [ustress,ustress_i] = max(rstress);
        ustrain = rstrain(ustress_i);
        %% Elastic Limit Approximation using Offset Yield Strength
        %When the stress-strain curve has no significant linear region, it
        %is recommended to use the offset yield strength to approximate the
        %linear elastic region of the material.
        %First Step is to apply a linear regression which pass through the
        %origin using the data from the first 20% strain (Trial and Error)
        offset = 2/100; %strain offset
        size_el = 20/100;
        [yield,fit_el,fit_off,el_mod] = offsetYield(rstrain,rstress,offset,size_el);
        %% Save Tensile Strength Data
        %Experimental Data
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).rstress = rstress;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).rstrain = rstrain;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).pstress = pstress;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).pstrain = pstrain;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).sstress = sstress;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).uStress = ustress; %Ultimate values
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).uStrain = ustrain;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yStress = yield(2); %Yield values
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yStrain = yield(1);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).elasticM = el_mod; %Elastic Modulus
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).elasticFit = fit_el;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).offsetFit = fit_off;
        
    end
end

end

function tensileStrength = getTensileStrengthBandsInd()
% Yellow, Red, Blue, Green, Black, Orange in acending order of thickness
load('TrainingSetProcStressV4.mat');
%Defining placeholders
f_names.mats =  fieldnames(PmatData); %contain the names of all the fields in PmatData
mats = length( f_names.mats ); %materials
tensileStrength = PmatData; %will store tensile data per test
for i = mats
    %Different tests were performed to each material. In here I find out
    %how many tests and their name.
    f_names.tests =  fieldnames(PmatData.(f_names.mats{i}));
    tests = length(  f_names.tests );
    for j = 1:tests
        %Inside the struct of each test there are many variables, as
        %follows: paths,pload,rload,rdis,rtime,stats and stress values
        f_names.data = fieldnames( PmatData.( f_names.mats{i} ).( f_names.tests{j} ));
        [~, P] = size( PmatData.( f_names.mats{i} ).( f_names.tests{j} ));
        %Clear previous data to store tensile properties in placeholder
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}) = [];
        for k=1:P
            rstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).rstress'; %rstress
            pstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).pstress'; %pstress
            sstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).sstress'; %sstress
            rstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).rstrain'; %rstrain
            pstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).pstrain'; %pstrain
            name = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).name'; %pstrain
            thickness = PmatData.( f_names.mats{i} ).(f_names.tests{j})(k).stats.thickness; %
            %Calculate Ultimate Values
            [ustress,ustress_i] = max(rstress);
            ustrain = rstrain(ustress_i);
            %% Elastic Limit Approximation using Offset Yield Strength            
            offset = 2/100; %strain offset
            size_el = 20/100;
            [yield,fit_el,fit_off,el_mod] = offsetYield(rstrain,rstress,offset,size_el);
            %% Save Tensile Strength Data
            %Experimental Data
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).rstress = rstress;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).rstrain = rstrain;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).pstress = pstress;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).pstrain = pstrain;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).sstress = sstress;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).uStress = ustress; %Ultimate values
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).uStrain = ustrain;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).yStress = yield(2); %Yield values
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).yStrain = yield(1);
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).elasticM = el_mod; %Elastic Modulus
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).elasticFit = fit_el;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).offsetFit = fit_off;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).name = name;
            tensileStrength.(f_names.mats{i}).(f_names.tests{j})(k).thickness = thickness;
        end
    end
end

end

function plotTensileStrength(tensileStrength)
%% Calculate size of structure
fontSize = 16;
f_names.N = fieldnames( tensileStrength);
N = length(f_names.N);
for i=1:N
    f_names.M = fieldnames( tensileStrength.(f_names.N{i}));
    M = length(f_names.M);
    f1 = figure;
    ax1 = gca;
    hold(ax1,'on');
    for j=1:M
        switch f_names.M{j}
            case "disR50"
                str_legends = "50 mm/min";
            case "disR250"
                str_legends = "250 mm/min";
            case "disR500"
                str_legends = "500 mm/min";
        end
        %Variables to plot
        x1 = tensileStrength.(f_names.N{i}).(f_names.M{j}).pstrain;
        y1 = tensileStrength.(f_names.N{i}).(f_names.M{j}).pstress;
        x2 = tensileStrength.(f_names.N{i}).(f_names.M{j}).rstrain;
        y2 = tensileStrength.(f_names.N{i}).(f_names.M{j}).rstress;
        p1 = plot(ax1, x1, y1,...
            'DisplayName',strcat(str_legends," Proc."));
        p2 = plot(ax1,x2  + 50, y2,...
            'DisplayName',strcat(str_legends," Raw"));
        ustrain(j) = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStrain;
        ustress(j) = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStress;
        p1.LineWidth = 2;
        p2.LineWidth = 2;
    end
    ax = gca;
    ax.FontSize = fontSize;
    grid on
    grid minor
    legend;
    xlabel('Displacement mm','FontSize',fontSize);
    ylabel('Load (N)','FontSize',fontSize);
    xlim( [ 0 max( ustrain ) ]);
    ylim( [ 0 max( ustress ) ])
    clear udis uload str_legends
    hold off
end
end

function plotOffsetYieldStrength(tensileStrength,f_names)
fontSize = 12;
LineWidth = 2;
f_names.N = fieldnames( tensileStrength);
N = length(f_names.N);

for i=1:N
    f_names.M = fieldnames( tensileStrength.(f_names.N{i}));
    M = length(f_names.M);   
    
    for j=1:M
        % %Chosing the title strings
        switch( f_names.M{j} )
            case {'disR50','disR50New'}
                f_names.tests_name = "50 mm/min";
            case {'disR250','disR250New'}
                f_names.tests_name = "250 mm/min";
            case {'disR500','disR500New'}
                f_names.tests_name = "500 mm/min";
        end
        %% Variables to plot
        pstrain = tensileStrength.(f_names.N{i}).(f_names.M{j}).pstrain;
        pstress = tensileStrength.(f_names.N{i}).(f_names.M{j}).pstress;
        sstress = tensileStrength.(f_names.N{i}).(f_names.M{j}).sstress;
        elasticFit = tensileStrength.(f_names.N{i}).(f_names.M{j}).elasticFit;
        offsetFit = tensileStrength.(f_names.N{i}).(f_names.M{j}).offsetFit;
        yStress = tensileStrength.(f_names.N{i}).(f_names.M{j}).yStress;
        yStrain = tensileStrength.(f_names.N{i}).(f_names.M{j}).yStrain;
        uStress = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStress;
        uStrain = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStrain;
        elasticM = tensileStrength.(f_names.N{i}).(f_names.M{j}).elasticM;
        
        f = figure();
        colors = get(gca,'ColorOrder');
        %Plots
        plot(pstrain*100,sstress*1e-6,'LineWidth',LineWidth,'Color',colors(1,:));
        hold on
        plot(pstrain(1:length(elasticFit))*100 , elasticFit*1e-6,'--',...
            'LineWidth',LineWidth,'Color',colors(2,:));
        plot(pstrain(1:length(offsetFit))*100 , offsetFit*1e-6,':',...
            'LineWidth',LineWidth,'Color',colors(3,:));
        plot(yStrain*100,yStress*1e-6,'go',...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',[1 1 1],...
            'MarkerSize',8,'LineWidth',LineWidth);
        hold off
        ylabel('Stress (MPa)','FontSize',fontSize);
        coor = find(pstrain*100<=100);
        ylim([0 sstress( coor(end) )*1e-6]);
        xlim([0 100]);
        grid on
        grid minor        
        ax = gca;
        ax.FontSize = fontSize;
        sTitle = strcat ( f_names.mats_name{i}," - ", f_names.tests_name);
        title(sTitle,'FontSize',fontSize);
        xlabel('Strain %','FontSize',fontSize);
        legend({strcat("Ultimate Strength [", num2str(round(uStress*1e-6,1))," MPa, ",num2str(round(uStrain*100))," %]"),...
            strcat("Modulus of Elasticity @ 20% [", num2str(round(elasticM*1e-6,1)), " MPa]"),...
            strcat("Offset Line @ 2% "),...
            strcat("Offset Yield Strength [",num2str(round(yStress*1e-6,1))," MPa, " ,num2str(round(yStrain*100))," %]")},...
            'Location','southeast','FontSize',fontSize);
        saveas(gcf,strcat( f_names.N{i} ,'_', f_names.M{j} ),'epsc');
    end
end



end

