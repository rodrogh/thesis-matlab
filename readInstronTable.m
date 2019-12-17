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
        for i=1:N
            %Ignore data after failure, i.e. after max load is found
            [~ , L] = max( dataStruct(i).load );
            load(i,1:L) = dataStruct(i).load(1:L);
            dis(i,1:L) = dataStruct(i).dis(1:L);
            time(i,1:L) = dataStruct(i).time(1:L);
            %Store length of each test
            test_lengths(i) = L;
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Statistic Values from the set of specimesn, this don't
        %%%%%% represent the ultimate displacement and load
        tempMax = max(load,[],2,'omitnan');
        minLoad = min(tempMax,[],1,'omitnan');
        maxLoad = max(tempMax,[],1,'omitnan');
        meanLoad = mean(tempMax,1,'omitnan');
        tempMax = max(dis,[],2,'omitnan');
        minDis = min(tempMax,[],1,'omitnan');
        maxDis = max(tempMax,[],1,'omitnan');
        meanDis = mean(tempMax,1,'omitnan');
        %Extra feature: Save the tes_lengths variable in statVals        
        statVals = struct('minLoad',minLoad,...
            'maxLoad',maxLoad,...
            'meanLoad',meanLoad,...
            'minDis',minDis,...
            'maxDis',maxDis,...
            'meanDis',meanDis);            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Save unified version of unprocessed data. The extraction of the
        %same amount of data for each test is required because tests were
        %sampled at different rates.
        [load,dis,time,test_lengths] = unifyData(N,load,dis,time,test_lengths);
        rload = mean(load,1); 
        rdis = mean(dis,1);
        rtime = mean(time,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Filtering slack data - slackFiltering(range,N,load,dis,time,plot_flag);
        [load,dis,time] = slackFiltering(100,N,load,dis,time,false);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Correcting negative offset found in data caused by calibrating the
        %load cell after the specimen was hold in place
        [load,dis,time] = negOffsetComp(N,load,dis,time);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The previous processes add padded zeros to the arrays, data must
        %be resized again. Then, unified by extracting the same amount of
        %data from all tests         
        [load,dis,time,test_lengths] = unifyData(N,load,dis,time,test_lengths);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Average all test into a single variable.
        %From here onwards, variables contain a single stress-strain curve
        pload = mean(load,1); %As long as dim is specified, this works with both methods
        pdis = mean(dis,1);
        ptime = mean(time,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Uncoment below section to apply smoothing to the data
        dataSize = length(pload); %updated size        
%         Find the suitable window size which causes a 5% attenuation of
%         the maximum value defined to observe
%         Define window size array
        window_t = unique(round( (linspace( 1, 10, 10)/100)*dataSize ));
        window_i =  find( window_t );
        window = window_t( window_i );
%         Test different window sizes from 1% of data to 10% in steps of 1%
%         Smoothing algorithm might benefit from padded data at the edges of
%         the curve
        window = round( ( linspace(1,10,10)*dataSize ) / 100);
        for n=1:length(window)
            sload(n,:) = smoothdata(pload,'sgolay',window(n) );
            %Measurement of how the maximum stress value is attenuated
%             temp_gof = diff(pload(n,:), 2 )./diff(rdis,2);
%             gof(n) = mean( temp_gof( ~isinf(temp_gof) ) );
            gof(n) = sqrt(mean((sload(n,:) - pload).^2 )); %REVISE THIS, tendency for upper boundary!!!!!!!!!!!!!!!!!!!!!!
        end       
        %Find desired window size
        [~ , index_best] = min( gof );        
        sload = sload(index_best,:);
        %Displaying useful data
%         disp('Savitzky-Golay');
%         disp(gof(index_best(end)));
%         disp('Window (%)');
%         disp(window(index_best)/length(pload) * 100);
%         disp('Window (size)');
%         disp(window(index_best));
        
%         figure
% %         subplot(1,2,1);
% %         yyaxis left
%         plot(rdis,rload,...
%             rdis,pload);
% %         yyaxis right
% %         plot(rdis(2:end),slope);
%         grid on;
% %         legend('Experimental','Savitzky-Golay');
%         title('Noise Removal');
%         xlim([0 0.1*max(rdis)]);
%         xlabel('Displacement (mm)');
%         ylabel('Load (N)');        
        
%         subplot(1,2,2);
%         plot(window/length(pload) * 100, gof,...
%             window/length(pload) * 100, (window./window)*mean(gof))
%         grid on;
%         title('Noise Removal');
%         xlabel('Window Size (%)');
%         ylabel('Mean Slope Variation (N/mm)');
%         %                 xlim([0 0.1*max(window)]);

    case 'Stress Relaxation'
        %Ensure deformation between current point and previous point is not
        %changing (Machine stable point)
        for i=1:N
            %Check first 10 data points, equilibrium must be reached by
            %then
            for j=1:r
                if dataStruct(i).dis(j)-dataStruct(i).dis(j+1) == 0
                    L = size(dataStruct(i).load(j:end));
                    load(i,1:L) = dataStruct(i).load(j:end);
                    dis(i,1:L) = dataStruct(i).dis(j:end);
                    time(i,1:L) = dataStruct(i).time(j:end);
                    %Statistical values of mean, min, max for ultimate load and
                    %displacement
                    maxLoad(i) = max(load(i,:));
                    maxDis(i) = max(dis(i,:));
                    break;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% Statistic Values
        minLoad = min(maxLoad);
        maxLoad = max(maxLoad);
        meanLoad = mean(maxLoad);
        minDis = min(maxDis);
        maxDis = max(maxDis);
        meanDis = mean(maxDis);
        statVals = struct('minLoad',minLoad,...
            'maxLoad',maxLoad,...
            'meanLoad',meanLoad,...
            'minDis',minDis,...
            'maxDis',maxDis,...
            'meanDis',meanDis);        
        mload = mean(load,1,'omitnan');
        mdis = mean(dis,1,'omitnan');
        mtime = mean(time,1,'omitnan');
        %Averaged raw data placeholders. Get rid of NaNs.
        nonNaN = find(isnan(mload)~=1);
        pdis = mdis(nonNaN);
        ptime =  mtime(nonNaN);
        pload =  mload(nonNaN);
        dataSize = length(pload);
        
        %         Find the suitable window size which causes a 5% attenuation of
        %         the maximum valuedefined to observe
        %         Define window size array
        window = linspace( 0.1, 20, 200);
        for n=1:length(window)
            pload(n,:) = smoothdata(pload,'sgolay',(window(n)/100)*dataSize);
            %             ploadg(n,:) = smoothdata(rload,'gaussian',(window(n)/100)*dataSize);
            %Measurement of how the maximum stress value is attenuated
            gof(n) =100 - 100*max(pload(n,:))./max(pload);
            %             gofg(n) =100 - 100*max(ploadg(n,:))./max(rload);
            
        end
        %Find desired window size
        index_best = find(gof <= 5);
        if isempty(index_best)
            index_best = 1;
        else    
            pload = pload(index_best(end),:);
        end
        disp(gof(index_best(end)));
        
%         figure
%         plot(rtime,rload,...
%             rtime,pload);
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
%Normalize the amount of data extracted from each test to have a
%uniform length along all tests.
%Resizing placeholders (Only required for method 0 but works for
%both methods)
load = NaN(N,min(test_lengths));
dis = NaN(N,min(test_lengths));
time = NaN(N,min(test_lengths));
for i=1:N
    %Linearly spaced vector of indeces based on size of smallest test
    temp4 = linspace(z_length(i) , test_lengths(i) , min(test_lengths) );
    temp4 = round(temp4);
    %Store unified data in placeholder
    load(i,:) = temp1( i,temp4  );
    dis(i,:) = temp2( i,temp4  );
    time(i,:) = temp3( i,temp4  );
    test_lengths(i) = length( load(i,:) );
end
uload = load;
udis = dis;
utime = time;
end

function [cload,cdis,ctime] = negOffsetComp(N,load,dis,time)
for i=1:N
    [neg_offset,neg_i] = min(  load(i,1:100), [],2 );
    if neg_offset < 0
        load(i,:) = load(i,:) + abs(neg_offset);
        load(i, 1:neg_i) = 0;
        dis(i, neg_i:end) = dis(i,1:end-(neg_i-1));
        dis(i, 1:neg_i) = 0;
        time(i, neg_i:end) = time(i,1:end-(neg_i-1));
        time(i, 1:neg_i) = 0;
    end
end
cload = load;
cdis = dis;
ctime = time;
end

function [cload,cdis,ctime] = slackFiltering(range,N,load,dis,time,plot_flag)
for i=1:N
    %Positive peaks with minimun height
    %TF are indices of peaks, represented by a 1
    %P are prominences of peaks
    [TF,P] = islocalmax(load(i,1:range));
    P = P(TF); %Use prominence of first peak as refrence to filter smaller peaks
    [TF_maxpeak,~] = islocalmax(load(i,1:range),'MinProminence',P(1));    
    maxpeak_i = find(TF_maxpeak > 0);
    maxpeak = load(i,maxpeak_i(1));
    %Finding negative peak
    [TF1,~] = islocalmin(load(i,1:range));
    minpeaks_i = find( TF1(1:end) == 1 );
    temp = NaN(1,range);
    temp(minpeaks_i) = load(i,minpeaks_i);
    temp(1:maxpeak_i(1)) = nan; %Delete negpeaks before first maxpeak
    [minpeak,minpeak_i] = min( temp ); %Find the smallest peak
    %Calculate decay percentage relative to maxpeak
    slack = 100*( abs(minpeak - maxpeak)/abs(maxpeak)  );
    uload = load(i,:);
    udis = dis(i,:);
    if abs(slack) > 30        
        %Resize data to filter slack section
        load(i,1:minpeak_i) = 0;
        dis(i,minpeak_i+1:end) = dis(i, 1:end-minpeak_i );
        dis(i,1:minpeak_i) = 0;
        time(i,minpeak_i+1:end) = time(i, 1:end-minpeak_i );
        time(i,1:minpeak_i) = 0;
    end
    %Plotting
    if plot_flag
        figure
        subplot(1,2,1);
        plot(udis,uload,...
            udis(maxpeak_i),uload(maxpeak_i),'*r',...
            udis(minpeak_i),uload(minpeak_i),'*g');
        stitle = strcat("Uncorrected Curve. Decay Detected ",num2str(slack),"%");
        title(stitle);
        xlim([0 dis(i,range)]);
        
        subplot(1,2,2);
        plot(dis(i,:),load(i,:),...
            dis(i,maxpeak_i),load(i,maxpeak_i),'*r',...
            dis(i,minpeak_i),load(i,minpeak_i),'*g');
        stitle = strcat("Corrected Curve. Decay Detected ",num2str(slack),"%");
        title(stitle);
        xlim([0 dis(i,range)]);        
    end
end
cload = load;
cdis = dis;
ctime = time;
end