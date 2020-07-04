%This script tests the generalization capabilities of the PL-SLS model,
%using the strain-dependent stiffness extracte from the file
%PiecewieL_SLS.m.
%In there, the strain rate of 500mm/min is used for all materials except
%the SR where 50mm/min is used.
%The extracted ke, under 10 different tolerance, hence number of strain segments
%is used fo obtain the prediction gof for the other available strain rates per material.
clear all
close all
clc
%% Plot visualization
% Change default axes fonts.
set(groot,'DefaultAxesFontName', 'Times New Roman')
set(groot,'DefaultAxesFontSize', 18)
% Change default text fonts.
set(groot,'DefaultTextFontname', 'Times New Roman')
set(groot,'DefaultTextFontSize', 18)
% Extra
set(groot, 'DefaultLegendLocation', 'best')
set(groot, 'defaultLineLineWidth', 2)
%% Initialization
load('PL.mat');
load('../Experimental work/TrainingSetProcTrainingSetProcStressV4.mat'); %Use this path for newly processed data
N = length( fieldnames(PmatData) );
fields.N = fieldnames(PmatData);
%Load Data
%Extract PL_ke_S
%Apply PL on other strain rates
%Calculate RMSE
%Plot 3D

for i=1:N
    M = length( fieldnames( PmatData.(fields.N{i}) ) );
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    for j=1:M
        if strcmp(fields.N{i},"SR" )
            if strcmp(fields.M{j},"disR50")
                continue;
            end
        else
            if strcmp(fields.M{j},"disR500")
                continue;
            end
        end
        %% Load parameters into placeholders
        force = PmatData.( fields.N{i} ).( fields.M{j} ). sstress;
        dis = PmatData.( fields.N{i} ).( fields.M{j} ). pstrain;
        time = PmatData.( fields.N{i} ).( fields.M{j} ). ptime;       
        dt = time(2) - time(1);
        switch fields.N{i}
            case "SR"
                modelFitField = "disR50";
            otherwise
                modelFitField = "disR500";
        end
        %These variables have multiple elements depending on the range of
        %tolerances used
        toleranceSize = length(modelFit.( fields.N{i} ).( modelFitField ));
        tempgofK = zeros(1,toleranceSize); %Known
        tempgofU = zeros(1,toleranceSize); %Unknown
        segments = zeros(1,toleranceSize);
        tol = zeros(1,toleranceSize);
        %FR Tolerance chosen is 5
        for t=1:toleranceSize
            tempgofK(t) = modelFit.( fields.N{i} ).( modelFitField )(t).gof;
%             tempgofK(t) = modelFit.( fields.N{i} ).( modelFitField )(t).gofNMAD;
            tol(t) = modelFit.( fields.N{i} ).( modelFitField )(t).tol;        
            segments(t) = modelFit.( fields.N{i} ).(modelFitField )(t).segments;
            PL_ke = modelFit.( fields.N{i} ).( modelFitField)(t).PL_ke;
            PL_ke_S = modelFit.( fields.N{i} ).( modelFitField )(t).PL_ke_S;
            PL_ke_dis = modelFit.( fields.N{i} ).( modelFitField )(t).PL_ke_dis;
            divisions = modelFit.( fields.N{i} ).( modelFitField )(t).disIndex;
            km = modelFit.( fields.N{i} ).( modelFitField )(t).km;
            thao = modelFit.( fields.N{i} ).( modelFitField )(t).thao;
            Eo = modelFit.( fields.N{i} ).( modelFitField )(t).Eo;
            %It is expected that the length of the fitte Ke is not equal to
            %the length of the test of other strain rate. When this
            %happens, compensate the length
            if length(dis) > length(PL_ke_S)
                deltaSize = length(dis) - length(PL_ke_S);
                paddEnd = PL_ke_S(end).*ones(1,deltaSize);
                PL_ke_S = [PL_ke_S , paddEnd];
            end
            %% PIECEWISE LINEARIZATION
            force_PL = zeros(1,length(dis));
            for a=2:length(dis)
                currentPL = PL_ke_S( PL_ke_dis(:) >= dis(a) );
                if isempty( currentPL )
                    currentPL = PL_ke_S(1);
                end
                force_PL(a) = ((currentPL(1) + km) * (dis(a) - dis(a-1)) + ...
                    currentPL(1)*dis(a)*dt/thao + force_PL(a-1)) / (1 + dt/thao);
            end
%             tempgofU(t) = 100*sqrt(mean(power(...
%                 ( 1/mean(force,'omitnan') ) .* (force - force_PL),2), 'omitnan' ));
            tempgofU(t) = 100*sqrt( mean( (force - force_PL).^2, 'omitnan')) ...
                ./ sqrt( mean(force.^2,'omitnan'));
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure
%             colors = get(gca,'ColorOrder');
%             plot(dis,force*1e-6,'Color',colors(6,:),'LineWidth',2);
%             hold on
%             plot(dis,force_PL*1e-6,':','Color',colors(7,:),'LineWidth',2);
%             hold off;
%             title( strcat( fields.N{i} ," - " ,  fields.M{j}) );
%             legend('Experimental Data','PL-SLS Fit');
%             xlabel('Strain');
%             ylabel('Stress (MPa)');
%             grid on;
%             ax = gca; % current axes
%             ax.Box = 'off';
%             ax.XMinorTick = 'on';
%             ax.YMinorTick = 'on';
%             ax.XMinorGrid = 'on';
%             ax.YMinorGrid = 'on';
%             ax.TickLength = [0.02 0.025];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        end
        
        gof.(fields.N{i}).(fields.M{j}) = tempgofU;
        gof.(fields.N{i}).(modelFitField) = tempgofK;
        gof.(fields.N{i}).segments = segments;
        gof.(fields.N{i}).tol = tol;
    end
end
%Extraction of extra parameters
for i=1:N
    M = length( fieldnames( PmatData.(fields.N{i}) ) );
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    dataHolder = zeros(M,toleranceSize);
    for j=1:M
        dataHolder(j,:) = gof.(fields.N{i}).(fields.M{j});  
    end
    %Extract mean and best values
    [meanbestGOF, I] = min( mean(dataHolder,1) );
    meanGOFTol = gof.(fields.N{i}).tol(I);
    meanGOTSegments = gof.(fields.N{i}).segments(I);
    
    [indbestGOF, I] = min( min(dataHolder,[],1) );
    indGOFTol = gof.(fields.N{i}).tol(I);
    indGOTSegments = gof.(fields.N{i}).segments(I);
    
    gofExtra.(fields.N{i}).best = meanbestGOF;
    gofExtra.(fields.N{i}).bestTol = meanGOFTol;
    gofExtra.(fields.N{i}).bestSegments = meanGOTSegments;
    gofExtra.(fields.N{i}).indbest = indbestGOF;
    gofExtra.(fields.N{i}).indbestTol = indGOFTol;
    gofExtra.(fields.N{i}).indbestSegments = indGOTSegments;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Loop for plotting

for i=1:N
    figure
    ax1 = gca;
    switch fields.N{i}
        case "NatR"
            plotName = "NatPolR";
        case "Nat100R"
            plotName = "NatR";
        otherwise
            plotName = fields.N{i};
    end
    M = length( fieldnames( PmatData.(fields.N{i}) ) );
    fields.M = fieldnames( PmatData.(fields.N{i}) );
    hold(ax1,'on');
    for j=1:M
        switch( fields.M{j} )
            case {'disR50','disR50New'}
                test_title = "50 mm/min";
            case {'disR250','disR250New'}
                test_title = "250 mm/min";
            case {'disR500','disR500New'}
                test_title = "500 mm/min";
        end
        %%%%
        
        switch fields.N{i}
            case "Nat100R"
                %%%%%
                switch fields.M{j}
                    case "disR50"
                        figure
                        ax2 = gca;
                        plot(ax2,tol,gof.(fields.N{i}).(fields.M{j}),'DisplayName', test_title, 'LineWidth',1.5);
                        title(ax2, strcat( "PL-SLS Performance - ", plotName),'Interpreter','latex' );
                        xlabel(ax2,'Tolerance (%)');
                        ylabel(ax2,'NMRSE');
                        legend(ax2,'Interpreter','latex');
                        ax2.Box = 'off';
                        ax2.XMinorTick = 'on';
                        ax2.YMinorTick = 'on';
                        ax2.XGrid = 'on';
                        ax2.YGrid = 'on';
                        ax2.TickLength = [0.02 0.025];
                        ax2.XLim = [10 100];
                        ax2.YLim = [550 600];
                    otherwise
                        plot(ax1,tol,gof.(fields.N{i}).(fields.M{j}),'DisplayName', test_title, 'LineWidth',1.5);
                end
                %%%%%
            otherwise
                plot(ax1,tol,gof.(fields.N{i}).(fields.M{j}),'DisplayName', test_title, 'LineWidth',1.5);
        end
        %%%%
        
        legend(ax1,'Interpreter','latex');
        title(ax1, strcat( "PL-SLS Performance - ", plotName),'Interpreter','latex' );
        xlabel(ax1,'Tolerance (%)');
        ylabel(ax1,'NMRSE');
        ax1.Box = 'off';
        ax1.XMinorTick = 'on';
        ax1.YMinorTick = 'on';
        ax1.XGrid = 'on';
        ax1.YGrid = 'on';
        ax1.TickLength = [0.02 0.025];
        ax1.XLim = [10 100];  
        
        
    end
    
    saveas(gcf,strcat( plotName ,'_Generalization'),'png');
    
end
        