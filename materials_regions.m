%% Description:
%Scripts to identify the toe and elastic regions of all the soft materials.
clear all
clc
close all
%% Initialization
%Plot visualization
% Change default axes fonts.
set(groot,'DefaultAxesFontName', 'Times New Roman')
set(groot,'DefaultAxesFontSize', 12)
% Change default text fonts.
set(groot,'DefaultTextFontname', 'Times New Roman')
set(groot,'DefaultTextFontSize', 12)
% Extra
set(groot, 'DefaultLegendLocation', 'best')
set(groot, 'defaultLineLineWidth', 2)
%Raw Dataset
dirname = "../Experimental Work/";
fnameT = strcat(dirname,"TrainingSetProcTrainingSetProcStressV4.mat");
% fnameT = strcat(dirname,"TrainingSetProcTrainingSetProcStressIndV4.mat");
fnameS = strcat(dirname,"TrainingSetProcSRelSetProcStressV4.mat");
%Strings for plotting
f_names.mats_name = {'Ethylene Polypropylene R';
    'Fluorocarbon R';
    'Natural R with Polyester';
    'Nitrile R';
    'Polyethylene R';
    'Silicone R';
    'Natural R'};
%Plotting values
f_names.smallRegion = { [1 2];
    [1 2.5];
    [0.6 3];
    [1 2];
    [0.6 0.25];
    [1 2];
    [2 1.5]};
% Elastic Toe Region by observation
f_names.toeRegion = { [0 0.2];
    [0 0.2 ];
    [0 0.04];
    [0 0.2];
    [0 0.2];
    [0 0.2];
    [0 0.4]};

%% Calculating ultimate values and yield offset
[tensileStrength] = getTensileStrength(fnameT,f_names);
% [tensileStrengthInd] = getTensileStrengthBandsInd(fname);
%% Ploting Stress-Strain Curves
% plotTensileStrength(tensileStrength,f_names);
%% Ploting Offset Yield Strength
% plotOffsetYieldStrength(tensileStrength,f_names);
%% Export TensileStrength Properties to a table
tensileTable = getFormattedTable(tensileStrength,10:17);
% tensileTableInd = getFormattedTableInd(tensileStrengthInd);
% [X,Y,Z]=getThicknessSurface(tensileStrengthInd);

%% Get Stress Relaxation
[stressRelax] = getStressRelax(fnameS,f_names);
sRelTable = getFormattedTable(stressRelax,8:11);
%% Plot Stress Relaxation
plotStressRelax(stressRelax,f_names)
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

function tbl= getFormattedTable(testData,range)
fields.N = fieldnames(testData);
N = length(fields.N);
T_data = table();
T_test = table();
for i=1:N
    fields.M = fieldnames(testData.(fields.N{i}));
    M = length(fields.M);
    for j=1:M
        T_test = [T_test; fields.N(i),fields.M(j)];
        tempT = struct2table(testData.(fields.N{i}).(fields.M{j}),'AsArray',true);
        T_data = [ T_data ; tempT(1,range) ];
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

function [yield_off, yfit_el, xfit_el, yoffset_fit, xoffset_fit , E] = offsetYield(x,y,offset,size_el)

%This function was re-writed on 11/Jan/2020
s_range = find( x <= size_el);
range_length = length( x( x <= size_el ) );
%Linear fit with offset
%Approach using polyfit
p = polyfit(x( x <= size_el ),y( x <= size_el ),1);
m = p(1);
%Fit linear regression to first portion
yfit_el = m.*x + p(2) ;
xfit_el = x;
%Fit linear reg to offset 
yoffset_fit = m.*(x - offset) + p(2);
xoffset_fit = x;
%Search is done in the short array (different index)
interIndex = find( (y - yoffset_fit) <= 0 );
yieldIndex = interIndex(1);
%Coordinates (x,y) of the offset yield strength
yield_off = [x(yieldIndex), y(yieldIndex)];
%Resize yfit_offset and xfit_offset to match where intercept is
yoffset_fit = yoffset_fit(1:yieldIndex );
xoffset_fit = xoffset_fit(1:yieldIndex );
xfit_el = xfit_el(1:yieldIndex);
yfit_el = yfit_el(1:yieldIndex);
%Output variables
E = m; %Modulus of elasticity
end

function tensileStrength = getTensileStrength(fname,f_names)
%% Initialization
%Loading the tensile strength processed data. The file: TrainingSetProc.mat
%contains the averaged data of all the specimens used n a single type of
%tensile strenght test.
load(fname);%Stress in Pa, and Strain NOT in %
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
        rstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstress'*1e-6; %rstress
        pstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstress'*1e-6; %pstress
        sstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).sstress'*1e-6; %sstress
        rstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstrain'; %rstrain
        pstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstrain'; %pstrain
        %Extract raw engineering stress (median)
        eng_stress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).eng_stress*1e-6;
        eng_strain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).eng_strain;
        %Calculate Ultimate Values
        [ustress,ustress_i] = max(sstress);
        ustrain = pstrain(ustress_i);
        %% Elastic Limit Approximation using Offset Yield Strength
        %When the stress-strain curve has no significant linear region, it
        %is recommended to use the offset yield strength to approximate the
        %linear elastic region of the material.
        %First Step is to apply a linear regression which pass through the
        %origin using the data from the first 20% strain (Trial and Error)
        offset = 20/100; %strain offset
%         size_el = 1/100;
        %Apply a linear regression at the middle of the toe range identified
        %visually
        size_el = (f_names.toeRegion{i}(2) - f_names.toeRegion{i}(1));
        [yield, yfit_el, xfit_el, yoffset_fit, xoffset_fit , E] = offsetYield(pstrain,sstress,offset,size_el);
        %% Calculate the Elastic Modulus of the second portion of the curve
        %i.e. from the offset yield strain to the end of the curve
        E2_Xrange = pstrain( pstrain >= yield(1) );
        E2_Yrange = sstress( pstrain >= yield(1) );
        %Resize range to focus on the middle section of the curve
        E2_Xrange = E2_Xrange(1:round(end/2));
        E2_Yrange = E2_Yrange(1:round(end/2));
        Ep = polyfit( E2_Xrange,E2_Yrange,1);
        E2 = Ep(1);
        E2_fit = Ep(1)*pstrain + Ep(2);
%         figure
%         plot(pstrain,sstress)     
       
        %% Save Tensile Strength Data
        %Experimental Data
        %Rounding constant, 2 for exporting table, 4 for plotting
        rc = 4;
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).rstress = round(rstress,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).rstrain = round(rstrain,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).pstress = round(pstress,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).pstrain = round(pstrain,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).sstress = round(sstress,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yelasticFit = round(yfit_el,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).xelasticFit = round(xfit_el,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yoffsetFit = round(yoffset_fit,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).xoffsetFit = round(xoffset_fit,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).eng_stress = round(eng_stress,rc); %10
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).eng_strain = round(eng_strain,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).uStress = round(ustress,rc); %Ultimate values
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).uStrain = round(ustrain,rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yStress = round(yield(2),rc); %Yield values
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).yStrain = round(yield(1),rc);
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).elasticM = round(E,rc); %Elastic Modulus
        tensileStrength.(f_names.mats{i}).(f_names.tests{j}).elasticM2 = round(E2,rc); %Elastic Modulus Large
        
        
    end
end

end

function tensileStrength = getTensileStrengthBandsInd(fname)
% Yellow, Red, Blue, Green, Black, Orange in acending order of thickness
load(fname);
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

function plotTensileStrength(tensileStrength,f_names)
%% Calculate size of structure
f_names.N = fieldnames( tensileStrength);
N = length(f_names.N);
for i=1:N
    f_names.M = fieldnames( tensileStrength.(f_names.N{i}));
    M = length(f_names.M);
    f1 = figure;
    ax1 = gca;
    %[0.1300 0.1100 0.7750 0.8150] 
    ax2 = axes('Position', [0.64 0.21 0.25 0.25 ] );
    hold(ax1,'on');
    hold(ax2,'on');
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
%         y1 = tensileStrength.(f_names.N{i}).(f_names.M{j}).pstress;
%         x2 = tensileStrength.(f_names.N{i}).(f_names.M{j}).rstrain;
%         y2 = tensileStrength.(f_names.N{i}).(f_names.M{j}).rstress;
        y3 = tensileStrength.(f_names.N{i}).(f_names.M{j}).sstress;
        %Plots for comparing Raw, Processed and Smoothed
%         p1 = plot(ax1, x1, y1,'--',...
%             'DisplayName',strcat(str_legends," Proc."));
%         p2 = plot(ax1,x2, y2,...
%             'DisplayName',strcat(str_legends," Raw"));
        p3 = plot(ax1,x1, y3,...
            'DisplayName',str_legends);
        p4 = plot(ax2,x1, y3);
%             'DisplayName',strcat(str_legends," Smoothed"));      
        ustrain(j) = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStrain;
        ustress(j) = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStress;
    end
    
    %Axes 1 Properties
    ax1.YGrid = 'on';
    ax1.XGrid = 'on';
    ax1.XMinorTick = 'on';
    ax1.YMinorTick = 'on';
    legend(ax1,'Location','northwest');
    title(ax1, f_names.mats_name{i} );
    xlabel(ax1,'Strain (mm/mm)');
    ylabel(ax1,'Stress (MPa)');
    xlim(ax1, [ 0 inf ]);
    ylim(ax1, [ 0 inf ]);
    cXTickDelta = ax1.XTick(2) - ax1.XTick(1);
    cYTickDelta = ax1.YTick(2) - ax1.YTick(1);
    ax1.XTick = [ax1.XTick (ax1.XTick(end)+cXTickDelta) ];
    ax1.YTick = [ax1.YTick (ax1.YTick(end)+cYTickDelta) ];  
    xlim(ax1, [ 0 ax1.XTick(end) ]);
    ylim(ax1, [ 0 ax1.YTick(end) ]);
    
%     xlim(ax1, [ 0 max( ustrain ) ]);
%     ylim(ax1, [ 0 max( ustress ) ]);

    %Axes 2 Properties
    ax2.YMinorGrid = 'on';
    ax2.XMinorGrid = 'on';
    ax1.XMinorTick = 'on';
    ax1.YMinorTick = 'on';
    xlabel(ax2,'Strain (mm/mm)');
    ylabel(ax2,'Stress (MPa)');
    xlim(ax2, [ 0 f_names.smallRegion{i}(1) ]);
    ylim(ax2, [ 0 f_names.smallRegion{i}(2) ]);
    title(ax2,"Initial Section");    
    
    %Final plot highlighted the zoomed area
    hold off
    saveas(gcf,strcat( f_names.N{i} ,"_StressStrain" ),'png');
    clear udis uload str_legends
    
    
end
end

function plotOffsetYieldStrength(tensileStrength,f_names)
fontSize = 12;
fontDelta = 4;
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
        yfit_el = tensileStrength.(f_names.N{i}).(f_names.M{j}).yelasticFit;
        xfit_el = tensileStrength.(f_names.N{i}).(f_names.M{j}).xelasticFit;
        yoffset_fit = tensileStrength.(f_names.N{i}).(f_names.M{j}).yoffsetFit;
        xoffset_fit = tensileStrength.(f_names.N{i}).(f_names.M{j}).xoffsetFit;
        yStress = tensileStrength.(f_names.N{i}).(f_names.M{j}).yStress;
        yStrain = tensileStrength.(f_names.N{i}).(f_names.M{j}).yStrain;
        uStress = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStress;
        uStrain = tensileStrength.(f_names.N{i}).(f_names.M{j}).uStrain;
        elasticM = tensileStrength.(f_names.N{i}).(f_names.M{j}).elasticM;
        
        f = figure();
        colors = get(gca,'ColorOrder');
        %Plots
        plot(pstrain,sstress,'LineWidth',LineWidth,'Color',colors(1,:));
        hold on
        plot(xfit_el , yfit_el,'--',...
            'LineWidth',LineWidth,'Color',colors(2,:));
        plot(xoffset_fit , yoffset_fit,':',...
            'LineWidth',LineWidth,'Color',colors(3,:));
        plot(yStrain,yStress,'go',...
            'MarkerEdgeColor',colors(3,:),'MarkerFaceColor',[1 1 1],...
            'MarkerSize',8,'LineWidth',LineWidth);
        hold off
%         coor = find(pstrain <= f_names.smallRegion{i}(1) );
        coor = find(pstrain <= 1 );
        ylim([0 sstress( coor(end) )]);
%         xlim([0 f_names.smallRegion{i}(1)]);
        xlim([0 1]);
        ax = gca;
        cXTickDelta = ax.XTick(2) - ax.XTick(1);
        cYTickDelta = ax.YTick(2) - ax.YTick(1);
        ax.XTick = [ax.XTick (ax.XTick(end)+cXTickDelta) ];
        ax.YTick = [ax.YTick (ax.YTick(end)+cYTickDelta) ];  
        xlim(ax, [ 0 ax.XTick(end) ]);
        ylim(ax, [ 0 ax.YTick(end) ]);
        grid on
        grid minor        
        
        ax.FontSize = fontSize+fontDelta;
        ax.XMinorTick = 'on';
        ax.Box = 'off';
        sTitle = strcat ( f_names.mats_name{i}," - ", f_names.tests_name);
        title(sTitle,'FontSize',fontSize+fontDelta);
        xlabel('Strain','FontSize',fontSize+fontDelta);
        ylabel('Stress (MPa)','FontSize',fontSize+fontDelta);
        legend({strcat("$[\sigma_u , \varepsilon_u$]=[", num2str(round(uStress,1))," MPa, ",num2str(round(uStrain,2))," ]"),...
            strcat("$E$ [", num2str(round(elasticM,1)), " MPa]"),...
            strcat("Offset Line @ 0.2 "),...
            strcat("$[\sigma_{y}, \varepsilon_{y}]$=[",num2str(round(yStress,1))," MPa, " ,num2str(round(yStrain,2))," ]")},...
            'Location','southeast','FontSize',fontSize+fontDelta/2,'Interpreter','latex');
        saveas(gcf,strcat( f_names.N{i} ,'_', f_names.M{j} ),'png');
    end
end



end

function sressRelax = getStressRelax(fname,f_names)
%% Initialization
%Loading the stress relaxation test processed data. The file: StressRelaxProc.mat
%contains the averaged data of all the specimens used n a single type of
%test
specimenLength = 33; %mm
load(fname);%Stress in Pa, and Strain NOT in %
%Defining placeholders
f_names.mats =  fieldnames(PmatData); %contain the names of all the fields in PmatData
mats = length( f_names.mats ); %materials
sressRelax = PmatData; %will store tensile data per test
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
        sressRelax.(f_names.mats{i}).(f_names.tests{j}) = [];
        %% Ultimate Values
        %The averaged values for rubber bands are stored in different
        %variables        
        rstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstress'*1e-6; %rstress
        pstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstress'*1e-6; %pstress
        sstress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).sstress'*1e-6; %sstress
        rstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rstrain'; %rstrain
        pstrain = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).pstrain'; %pstrain
        rtime = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).rtime'; %pstrain
        ptime = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).ptime'; %pstrain
        %Extract initial stress stress (median)
        eng_stress = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).eng_stress*1e-6;
        %Extract stress at end of test
        eng_stress_end = PmatData.( f_names.mats{i} ).(f_names.tests{j})(1).eng_stress_end*1e-6;
        delta_stress = (eng_stress - eng_stress_end)/eng_stress;
        %% Save Stress Relaxation Data
        %Experimental Data
        %Rounding constant, 2 for exporting table, 4 for plotting
        rc = 4;
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).rstress = round(rstress,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).rstrain = round(rstrain,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).pstress = round(pstress,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).pstrain = round(pstrain,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).sstress = round(sstress,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).rtime = round(rtime,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).ptime = round(ptime,rc); %7
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).e_init = round(pstrain(1)*specimenLength,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).eng_stress = round(eng_stress,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).eng_stress_end = round(eng_stress_end,rc);
        sressRelax.(f_names.mats{i}).(f_names.tests{j}).delta_stress = round(delta_stress,rc);
    end
end

end

function plotStressRelax(stressRelax,f_names)
%% Calculate size of structure
f_names.N = fieldnames( stressRelax);
N = length(f_names.N);
fig1_180 = figure;
fig1_15 = figure;
fig2_180 = figure;
fig2_15 = figure;
for i=1:N
    f_names.M = fieldnames( stressRelax.(f_names.N{i}));
    M = length(f_names.M);
    for j=1:M
        %Variables to plot        
        x = stressRelax.(f_names.N{i}).(f_names.M{j}).ptime;
        y = stressRelax.(f_names.N{i}).(f_names.M{j}).sstress;
        %Plotting variables
        epsilon = round ( stressRelax.(f_names.N{i}).(f_names.M{j}).e_init);
        str_legends = strcat( "$\varepsilon_o$=",num2str( round(epsilon/33,2)), " ", f_names.mats_name{i} );
        switch f_names.M{j}            
            case "L5min"
                continue;
            case "L15min"                
                strTime = "15 min";
                if i == 5 || i == 7
                    figure(fig1_15)
                    if i==7
                        yyaxis right                        
                    else
                        yyaxis left
                    end
                else
                    figure(fig2_15)
                end
            case "L180min"
                strTime = "180 min";
                if i == 5
                    figure(fig1_180)                    
                else
                    figure(fig2_180)                                     
                end
        end
        ax1 = gca;
        hold(ax1,'on');
        plot(ax1,x/60,y,'DisplayName',str_legends); %Convert seconds to minutes
        %Axes 1 Properties
        ax1.YGrid = 'on';
        ax1.XGrid = 'on';
        ax1.XMinorTick = 'on';
        ax1.YMinorTick = 'on';
        legend(ax1,'Location','northeast','Interpreter','latex');
        strTitle = strcat( "Stress Relaxation ", strTime);
        title(ax1, strTitle,'Interpreter','latex');
        xlabel(ax1,'Time (minutes)');
        ylabel(ax1,'Stress (MPa)');
        hold off
        saveas(gcf,strcat( f_names.N{i} ,"_StressRel" ),'png');
        clear str_legends
    end
end
end