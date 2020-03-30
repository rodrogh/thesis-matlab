%Author: Rodrigo Daniel Solis Ortega
%The following function extracts the data from the .csv files created by
%the Instron Testing Machine.
%The function takes the following inputs:
% N -> Integer, defines how many .csv files to read
% iRow -> Integer, row coordinate where tha data inside the table begins
% iCol -> Integer, Column coordinate where tha data inside the table begins
% path -> String Cell Array, contains the paths where the files are located
%exp -> defines which processing approach to implement. 'TensileStrength' or
%'StressRelaxation'
%The data is processed here. The function outputs are:
%disMean -> Averaged displacement data of the N experiments
%loadMean -> Averaged load data of the N experiments
%time -> timestamp of the experiment
function [rload, rdis, rtime, pload, pdis, ptime, sload, statVals] ...
    = readInstronTable(N,iRow,iCol,paths,test)
%Define data holders
% rawData = NaN(N,r,c);
dataStruct = struct('load',{},'dis',{},'time',{});
%Variable for holding the length of each file
a = zeros(1,N);
N = length(paths);
%Read data from data tables
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [row,~] = size(temp);
    a(i) = row; %Save the length of each file (rows)
    %Save all rawdata
    %rawData(i,1:row,:) = temp;
    %Save desired raw data
    dataStruct(i).load = temp(:,3);
    dataStruct(i).dis = temp(:,2);
    dataStruct(i).time = temp(:,1);
    clear p;
    clear temp;
%     figure
%     plot(dataStruct(i).dis,dataStruct(i).load);
%     grid on;
%     title('Raw Data');
%     xlabel('Displacement (mm)');
%     ylabel('Load (N)');
end
%Concatenate data inside the struct into a variable of same length, using
%vector a. Some files will have less rows than others, these empty rows are
%filled with NaN
[r,rIndex] = max(a);
%For tensile strength, get the index in which the largest time vector is
%found
load = NaN(N,r);
dis = NaN(N,r);
time = NaN(N,r);
test_lengths = NaN(N,1);
switch (test)
    case 'Tensile Strength'
        %Algorithm Features
        %1. Filter data beyond the failure point.     
        %2. Extract the mean, min and max of the group of 5 set of data
        %3. Filter undesired slack stress data at the beginning of the
        %curve
        %4. Compensate negative offset
        %5. Filter padded zeros, added from previous steps and Extract 
        %an equal number of data point from all test prior to average them 
        %into a single dataset (method 0)
        %6. Calculate the mean value along all the data set to create a
        %single stress-strain curve 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Processing stage, incorporate the N experiments into one dataset
%         f_temp = figure;
%         ax_temp = gca;
%         ax_small = axes('Position',[0.2 0.65 0.25 0.25 ]);
%         hold(ax_temp,'on');
%         hold(ax_small,'on');
        for i=1:N
            %Plot Raw data
%             plot(ax_temp,dataStruct(i).dis,dataStruct(i).load);
%             plot(ax_small,dataStruct(i).dis,dataStruct(i).load);
            %Ignore data after failure, i.e. after max load is found
            [~ , L] = max( dataStruct(i).load );
            load(i,1:L) = dataStruct(i).load(1:L);
            dis(i,1:L) = dataStruct(i).dis(1:L);
            time(i,1:L) = dataStruct(i).time(1:L);
            %Store length of each test
            test_lengths(i) = L;
        end
%         ax_temp.YGrid = 'on';
%         ax_temp.XGrid = 'on';
%         ax_small.YMinorGrid = 'on';
%         ax_small.XMinorGrid = 'on';
%         ax_small.XMinorTick = 'on';
%         ax_small.YMinorTick = 'on';
%         xlabel(ax_temp,'Displacement (mm)');
%         ylabel(ax_temp,'Load (N)');
%         xlabel(ax_small,'Displacement (mm)');
%         ylabel(ax_small,'Load (N)');
%         title(ax_temp,'Tensile Strength Test Raw Data');
%         xlim(ax_small,[0 2]);
%         ylim(ax_small,[-1 1]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find displacement rate
        disR = max ( dis( time <= 1) );
        %Resample and redefine largest size
        %Find lowest sampling rate
        delta_time = time(:,2) - time(:,1);
        [max_time, max_time_i] = max(delta_time); %Find lowest sampling, i.e. highest period
        longest_time = max(time,[],2); %Most complete time vector
        %Define new sampling time vector
        sample_time = 0: max_time : max(longest_time);
        for i=1:N
            k=1; %New counter
            M = length(time(i,:));
            for j=1:M
               if time(i,j) == 0 || mod(time(i,j),max_time)==0
                    newLoad(i,k) = load(i,j);
                    newDis(i,k) = dis(i,j);
                    newTime(i,k) = time(i,j);
                    k=k+1;
               end
               %New test_lengths
               test_lengths(i) = k-1;
            end
        end
        %New data placeholders
        [~,new_i] = max(newTime,[],2);
        new_largest = max( new_i );
        load = NaN(N,new_largest);
        dis = NaN(N,new_largest);
        time = NaN(N,new_largest);
        for i=1:N
            load(i,1:new_i(i)) = newLoad(i,1:new_i(i));
            dis(i,1:new_i(i)) = newDis(i,1:new_i(i));
            time(i,1:new_i(i)) = newTime(i,1:new_i(i));            
        end
        %Make initial displacement equals zero to avoid unnecessary
        %processing
        dis(:,1) = 0;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Statistic Values from the set of specimesn, this don't
        %%%%%% represent the ultimate displacement and load
        tempMax = max(load,[],2,'omitnan');
        minLoad = min(tempMax,[],1,'omitnan');
        maxLoad = max(tempMax,[],1,'omitnan');
        meanLoad = mean(tempMax,1,'omitnan');
        medLoad = median(tempMax,1,'omitnan');
        tempMax = max(dis,[],2,'omitnan');
        minDis = min(tempMax,[],1,'omitnan');
        maxDis = max(tempMax,[],1,'omitnan');
        meanDis = mean(tempMax,1,'omitnan');
        medDis = median(tempMax,1,'omitnan');
        %Extra feature: Save the test_lengths variable in statVals        
        statVals = struct('minLoad',minLoad,...
            'maxLoad',maxLoad,...
            'meanLoad',meanLoad,...
            'medLoad',medLoad,...
            'minDis',minDis,...
            'maxDis',maxDis,...
            'meanDis',meanDis,...
            'medDis',medDis);            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Save unified version of unprocessed data. The extraction of the
        %same amount of data for each test is required because tests were
        %sampled at different rates. Discard uneven failure sections
%         [load,dis,time,test_lengths] = unifyData(N,load,dis,time,test_lengths);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Resize variables to discard the uneven failure points effect. Find
        %the smallest test without NaN values, this is stored in
        %test_lengths
        load = load(:,1:min(test_lengths));
        dis = dis(:,1:min(test_lengths));
        time = time(:,1:min(test_lengths));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Unify Raw data
        rload = mean(load,1); 
        rdis = mean(dis,1);
        rtime = mean(time,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Filtering slack data - slackFiltering(range,N,load,dis,time,plot_flag);
        initSection = find(rtime <= 2);
        initSectionSize = length(initSection);
        [load,dis,time] = slackFiltering(initSection,N,load,dis,time,false);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Correcting negative offset found in data caused by calibrating the
        %load cell after the specimen was hold in place
        [load,dis,time,test_lengths] = negOffsetComp(N,load,dis,time,test_lengths);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The previous processes add padded zeros to the arrays, data must
        %be resized again. Then, unified by extracting the same amount of
        %data from all tests         
%         [load,dis,time,test_lengths] = unifyData(N,load,dis,time,test_lengths);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Average all test into a single variable.
        %From here onwards, variables contain a single stress-strain curve
        pload = mean(load,1); %As long as dim is specified, this works with both methods
        pdis = mean(dis,1); %omitnan caused the abrupt changes at the end of the curve
        ptime = mean(time,1);
        %Find nans and resize
        nanVar = isnan(pload);
        pload = pload( ~nanVar );
        pdis = pdis( ~nanVar );
        ptime = ptime( ~nanVar );        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Correcting negative offset found in unified data caused by the 
        %mean calculation
        [pload,pdis,ptime,test_lengths] = negOffsetComp(1,pload,pdis,ptime,test_lengths);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Find nans and resize
        nanVar = isnan(pload);
        pload = pload( ~nanVar );
        pdis = pdis( ~nanVar );
        ptime = ptime( ~nanVar );           
        %Mean value of initial section
        initSection = find(ptime <= 1);
        initSectionSize = length(initSection); %Window for smoothing algorithm
        %Find peaks and calculate mean values. If different than 0. Then
        %smooth
        TFmax = nan(1,initSectionSize);
        TFmin = nan(1,initSectionSize);
        TFmax(:) = islocalmax(pload (initSection));
        TFmin(:) = islocalmin(pload (initSection));
        peakTopeak = mean( pload(TFmax==1)) - mean( pload(TFmin==1) );
        initSectionMean = mean( peakTopeak );
        
        if initSectionMean ~= 0            
            sload = smoothdata(pload,'sgolay', initSectionSize);           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Correcting negative offset found in unified data caused by the
            %smoothing            
            if sload(1:20)<0
               disp('------------------------NEGATIVE------------------------') 
            end
            [sload,pdis,ptime,test_lengths] = negOffsetComp(1,sload,pdis,ptime,test_lengths);
            if sload(1:20)<0
               disp('------------------------STILL NEGATIVE------------------------') 
            end
            %                 figure
            %                 plot(pdis,pload,'-','LineWidth',1);
            %                 hold on
            %                 plot(pdis,sload(n,:),'--','LineWidth',2);
            %                 hold off
            %                 grid on;
            %                 legend('Experimental','Savitzky-Golay');
            %                 title('Noise Removal');
            %                 xlim([0 0.1*max(rdis)]);
            %                 xlabel('Displacement (mm)');
            %                 ylabel('Load (N)');
            
%             figure
%             %         subplot(1,2,1);
%             %         yyaxis left
%             plot(rdis,rload,'-','LineWidth',1);
%             hold on
%             plot(pdis,pload,'--','LineWidth',1);
%             plot(pdis,sload,':');
%             hold off
%             %         yyaxis right
%             %         plot(rdis(2:end),slope);
%             grid on;
%             legend('Experimental','Processed','Savitzky-Golay');
%             title('Noise Removal');
%             xlim([0 0.1*max(rdis)]);
%             xlabel('Displacement (mm)');
%             ylabel('Load (N)');            
        else
            sload = pload;
            disp('Smoothing not required');
        end
        %Extra step to make sure that the rstrain and pstrain variables do
        % not have extra zeros at the beginning. disR is in mm
%         rdis_zeros = length( find( rdis == 0 ) );
%         if rdis_zeros > 1 && ~isempty(rdis_zeros)
%             disp("Zero corrected for rdis");
%             dis_step = rdis(end) - rdis(end-1);
%             rdis = circshift( rdis,-(rdis_zeros-1) );
%             for z=(rdis_zeros-1):-1:0
%                 rdis(end-z) = rdis((end-1)-z) + dis_step;            
%             end
%         end
%         pdis_zeros = length( find( pdis == 0 ) );
%         if pdis_zeros > 1 && ~isempty(pdis_zeros)
%             dis_step = rdis(end) - pdis(end-1);
%             pdis = circshift( pdis,-(pdis_zeros-1) );
%             pdis(end) = pdis(end-1) + dis_step;
%             disp("Zero corrected for pdis");
%             for z=(pdis_zeros-1):-1:0
%                 pdis(end-z) = pdis((end-1)-z) + dis_step;            
%             end
%         end
% rdis = disR*rtime;     
% pdis = disR*ptime;
%      
    case 'Stress Relaxation'
        %Resample the data with 1Hz frequency. Datasets are too long
        %Resample and redefine largest size
        %Find lowest sampling rate
        delta_time = 1;
        size_time = length( dataStruct(1).time);
        clear load dis time
        for i=1:N
            samplingFreq(i) = 1/(dataStruct(i).time(2) - dataStruct(i).time(1));
            k=0; %New counter
            for j=1:size_time-1
                if ( dataStruct(i).dis(j) - dataStruct(i).dis(j+1) ) == 0
                    if (mod( round(dataStruct(i).time(j) ),delta_time)==0 )
                        k=k+1;
                        newload(i,k) = dataStruct(i).load(j);
                        newdis(i,k) = dataStruct(i).dis(j);
                        newtime(i,k) = dataStruct(i).time(j) ;
                    end
                end
               %New test_lengths
               test_lengths(i) = k;
            end
            %Due to large rounding for scanning properly, unique must be used
            [~,uniqueTime(i,:),~] = unique( round ( newtime(i,:) ) );
            load(i,:) = newload(i,uniqueTime(i,:) );
            dis(i,:) = newdis(i,uniqueTime(i,:) );
            time(i,:) = newtime(i,uniqueTime(i,:) );
        end
        
        %New data placeholders
%         [~,new_i] = max(newTime,[],2);
%         new_largest = max( new_i );
%         load = NaN(N,new_largest);
%         dis = NaN(N,new_largest);
%         time = NaN(N,new_largest);
        %Ensure deformation between current point and previous point is not
        %changing (Machine stable point)
%         for i=1:N
%             %Check first 10 data points, equilibrium must be reached by
%             %then
%             for j=1:r
%                 if dataStruct(i).dis(j)-dataStruct(i).dis(j+1) == 0
%                     L = size(dataStruct(i).load(j:end));
%                     load(i,1:L) = dataStruct(i).load(j:end);
%                     dis(i,1:L) = dataStruct(i).dis(j:end);
%                     time(i,1:L) = dataStruct(i).time(j:end);
%                     %Statistical values of mean, min, max for ultimate load and
%                     %displacement
%                     maxLoad(i) = max(load(i,:));
%                     maxDis(i) = max(dis(i,:));
%                     break;
%                 end
%             end
%         end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Statistic Values
        minLoad = median(  min ( load,[],2 ) ); %Load at the end of each test
        maxLoad = max( max ( load,[],2 ) );
        meanLoad = mean( max ( load,[],2 ) );
        medLoad = median( max ( load,[],2 ) );
        minDis = min( max ( dis,[],2 ) );
        maxDis = max( max ( dis,[],2 ) );
        meanDis = mean( max ( dis,[],2 ) );
        medDis = median( max ( dis,[],2 ) );
        statVals = struct('minLoad',minLoad,...
            'maxLoad',maxLoad,...
            'meanLoad',meanLoad,...
            'medLoad',medLoad,...
            'minDis',minDis,...
            'maxDis',maxDis,...
            'meanDis',meanDis,...
            'medDis',medDis);        
        rload = mean(load,1,'omitnan');
        rdis = mean(dis,1,'omitnan');
        rtime = mean(time,1,'omitnan');
        pload = rload;
        pdis = rdis;
        ptime = rtime;
        %There is no nan data on placeholders due to sampling performed
        %above
        dataSize = length(pload);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Mean value of initial section
        initSection = find(ptime <= 50);
        initSectionSize = length(initSection); %Window for smoothing algorithm
        %Find peaks and calculate mean values. If different than 0. Then
        %smooth
        TFmax = nan(1,initSectionSize);
        TFmin = nan(1,initSectionSize);
        TFmax(:) = islocalmax(pload (initSection));
        TFmin(:) = islocalmin(pload (initSection));
        peakTopeak = mean( pload(TFmax==1)) - mean( pload(TFmin==1) );
        initSectionMean = mean( peakTopeak );
        if initSectionMean ~= 0
            sload = smoothdata(pload,'sgolay', 50 );
%             %             ploadg(n,:) = smoothdata(rload,'gaussian',(window(n)/100)*dataSize);
%             %Measurement of how the maximum stress value is attenuated
%             gof(n) =100 - 100*max(sload(n,:))./max(pload);
%             %             gofg(n) =100 - 100*max(ploadg(n,:))./max(rload);
        else
            sload = smoothdata(pload,'sgolay', 50 );
            sload = pload;
            disp("Smoothing not required");
        end
%         %Find desired window size
%         index_best = find(gof <= 5);
%         if isempty(index_best)
%             index_best = 1;
%         else    
%             sload = sload(index_best(end),:);
%         end
%         disp(gof(index_best(end)));
%         
%         figure
%         plot(ptime,pload,...
%             ptime,sload);
%         grid on;
%         %                 xlim([0 1000]);
%         legend('Experimental','Savitzky-Golay');
%         title('Noise Removal');
%         xlabel('Time (seconds)');
%         ylabel('Load (N)');
        
        %Plot fot comparison analysis between two methods
        %         figure
        %         plot(window,gof, window,gofg);
        %         grid on;
        %         title('Savitzky-Golay vs Gaussian-weighted');
        %         xlabel('Windows Size (% of total data)');
        %         ylabel('Attenuation (%)');
        %         legend('Savitzky-Golay','Gaussian-weighted');
        
    otherwise
        disp('Unkown method');
end
end

function [uload,udis,utime,test_lengths] = unifyData(N,load,dis,time,test_lengths)
%This function also updates the length of each test. This is useful is 
%the function is called prior to negOffsetComp        
temp1 = load;
temp2 = dis;
temp3 = time;
z_length = ones(N,1);
%Find zeros at the beginning of data and update useable length
for i=1:N
    z = find( temp1(i,:) == 0 );
    if ~isempty(z)
        z_length(i) = length(z);
        test_lengths(i) = test_lengths(i) - (z_length(i)-1);
    end
end
%% Discard failure section with different endings. Take the smallest of all
% the test as the reference.
%Normalize the amount of data extracted from each test to have a
%uniform length along all tests. 
%Resizing placeholders (Only required for method 0 but works for
%both methods)
% load = NaN(N,min(test_lengths));
% dis = NaN(N,min(test_lengths));
% time = NaN(N,min(test_lengths));
load = temp1(:,1:min(test_lengths));
dis = temp2(:,1:min(test_lengths));
time = temp3(:,1:min(test_lengths));
% for i=1:N    
%     %Linearly spaced vector of indeces based on size of smallest test
%     temp4 = linspace(z_length(i) , test_lengths(i) , min(test_lengths) );
%     temp4 = round(temp4);
%     %Store unified data in placeholder
%     load(i,:) = temp1( i,temp4  );
%     dis(i,:) = temp2( i,temp4  );
%     time(i,:) = temp3( i,temp4  );
%     test_lengths(i) = length( load(i,:) );
% end
uload = load;
udis = dis;
utime = time;
end

function [cload,cdis,ctime,test_lengths] = negOffsetComp(N,load,dis,time,test_lengths)
neg_offset = zeros(N,1);
neg_i = zeros(N,1);
for i=1:N
    %Search all tests for negative offsets
    [neg_offset(i),neg_i(i)] = min(  load(i,:), [],2 );
    %The largest neg_i found will dictates the new size for the unified
    %test arrays
    %It could be the case that the test does not have a negOffset
end
%Make sure that at least one of the minimum found is negative, otherwise,
%exit function
if min(neg_offset) < 0
    %Define new placeholders
    maxNeg_i = max(neg_i);
    oldSize = length( load(1,:) );
    newSize = oldSize - maxNeg_i;
    newload = load(:,1:newSize);
    newdis = dis(:,1:newSize);
    newtime = time(:,1:newSize);
    test_lengths(i) = newSize;
    for i=1:N
        %Compensate negOffset if negative
        if neg_offset(i) < 0
            load(i,:) = load(i,:) + abs(neg_offset(i) );
            %Make the point in which the negOffset was found, the new starting
            %point of the curve
            load(i,:) = circshift( load(i,:), -neg_i(i) );
            newload(i,:) = load(i,1:newSize);
            %The variables for time and dis do no need to be shifted. In this
            %way, they both will represent the beginning of the test
            newdis(i,:) = dis(i,1:newSize);
            newtime(i,:) = time(i,1:newSize);        
        end
        
        %     while neg_offset < 0
        %         %Compensate negOffset
        %         load(i,:) = load(i,:) + abs(neg_offset);
        %         dis(i, neg_i:end) = dis(i,1:end-(neg_i-1));
        %         dis(i, 1:neg_i) = 0;
        %         time(i, neg_i:end) = time(i,1:end-(neg_i-1));
        %         time(i, 1:neg_i) = 0;
        %         [neg_offset,neg_i] = min(  load(i,:), [],2 );
        %     end
    end
    cload = newload;
    cdis = newdis;
    ctime = newtime;
else
    cload = load;
    cdis = dis;
    ctime = time;
end
end



function [cload,cdis,ctime] = slackFiltering(range,N,load,dis,time,plot_flag)
rangeSize = length(range);
for i=1:N    
    %Positive peaks with minimun height
    %TF are indices of peaks, represented by a 1
    %P are prominences of peaks
    [TF,P] = islocalmax(load(i,range));
    if TF==0
       disp('No take-up slack detected');
       cload = load;
       cdis = dis;
       ctime = time;
       return; 
    end
    P = P(TF); %Use prominence of first peak as refrence to filter smaller peaks
    [TF_maxpeak,~] = islocalmax(load(i,range),'MinProminence',P(1));    
    maxpeak_i = find(TF_maxpeak > 0);
    maxpeak = load(i,maxpeak_i(1));
    %Finding negative peak
    [TF1,~] = islocalmin(load(i,range));
    minpeaks_i = find( TF1(1:end) == 1 );
    temp = NaN(1,rangeSize);
    temp(minpeaks_i) = load(i,minpeaks_i);
    temp(1:maxpeak_i(1)) = nan; %Delete negpeaks before first maxpeak
    [minpeak,minpeak_i] = min( temp ); %Find the smallest peak
    %Calculate decay percentage relative to maxpeak
%     slack = 100*( abs(minpeak - maxpeak)/abs(maxpeak)  );
    slack = 100*( (maxpeak - minpeak)/maxpeak  );
    uload = load(i,:);
    udis = dis(i,:);
    if abs(slack) > 30        
        %Resize data to filter slack section
        load(i,:) = circshift(load(i,:),-minpeak_i);
        load(i,end-minpeak_i:end) = nan;
%         dis(i,minpeak_i+1:end) = dis(i, 1:end-minpeak_i );
%         dis(i,1:minpeak_i) = 0;
%         time(i,minpeak_i+1:end) = time(i, 1:end-minpeak_i );
%         time(i,1:minpeak_i) = 0;
    else
        disp('No take-up slack detected');
        cload = load;
        cdis = dis;
        ctime = time;
        return;
    end
    %Plotting
%     plot_flag=true;
    if plot_flag
        figure
        subplot(1,2,1);
        plot(udis,uload,...
            udis(maxpeak_i),uload(maxpeak_i),'*r',...
            udis(minpeak_i),uload(minpeak_i),'*g');
        stitle = strcat("Uncorrected Curve. Decay Detected ",num2str(slack),"%");
        title(stitle);
        xlim([0 dis(i,range(end))]);
        
        subplot(1,2,2);
        plot(dis(i,:),load(i,:),...
            dis(i,maxpeak_i),load(i,maxpeak_i),'*r',...
            dis(i,minpeak_i),load(i,minpeak_i),'*g');
        stitle = strcat("Corrected Curve. Decay Detected ",num2str(slack),"%");
        title(stitle);
        xlim([0 dis(i,range(end))]);        
    end
end
cload = load;
cdis = dis;
ctime = time;
end