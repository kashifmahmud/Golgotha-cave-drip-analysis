clear; home
clc;
addpath('/Users/kashifmahmud/UNSW/Golgotha drip data')  

%% AWAP data import fro HEP
[AWAPm,Dm] = xlsread('awap2-operational-evaporation-timeseries-monthly-Golgotha with Andys edits.xlsx'); % Import HEP calculations
[AWAPw,Dw] = xlsread('awap2-operational-evaporation-timeseries-weekly-Golgotha.xlsx'); % Import HEP calculations
% [AWAPd Dd] = csvread('awap-daily-rainfall-golgotha.csv'); % Import daily rainfall data
[AWAPd,Dd] = xlsread('awap-daily-rainfall-golgotha.xlsx'); % Import daily rainfall data
rainfalld = AWAPd(:,3);
% rainfalld = rainfalld(2984:end,1);
rainfall = AWAPw(:,15);
rainfall = rainfall(51:end,1);
P_ETw = AWAPw(:,17);
P_ETw = P_ETw(55:end,1);
P_ETm = AWAPm(:,27);
P_ETm = P_ETm(end-35:end,1);

% Recharge = AWAPm(:,22);
% TF = isnan(Recharge);
% Recharge(TF==1) = [];
% % Recharge(:) = Recharge(:)/min(Recharge);
% Recharge_norm = (Recharge - min(Recharge)) / ( max(Recharge) - min(Recharge) );
% Dm(:,2:end) = [];
% Recharge_time = Dm(TF==0);
% formatIn = 'dd/mm/yyyy';
% Recharge_time = datenum(Recharge_time(:),formatIn);

%% Rainfall calculation
startDate = datenum('05-Apr-2012');
endDate = datenum('06-Mar-2015');
% increment = datenum('00-01-0000');
rainfall(:,2) = linspace(startDate,endDate,size(rainfall,1));

%%
[AWAP_avg D_avg] = xlsread('Rainfall_monthly_avg_data.xlsx'); % Import monthly rainfall data
% rain_avg = AWAP_avg(18:end, [5:13 2:13 2:13 2:4]);
rain_avg = AWAP_avg(18:end, [5:13 2:4]); 

formatIn_m = 'mm/yyyy';
startDate_m = datenum('Apr-2012');
endDate_m = datenum('Mar-2014');
rainfallm_2012(:,1) = [80.2	127.8	265.2	86.2	144.8	170.4	31.8	61.6	23.0 10.0	7.6	49.4]';
rainfallm_2013(:,1) = [23.0	220.8	189.0	251.2	184.2	218.4	48.4	33.6	4.2 4.4	1.4	29.6]';
rainfallm_2014(:,1) = [24.0	149.6	180.6	213.6	168.4	88.6	32.6	50.0	1.0    0.0	 13.8	56.2]';
rainfallm(:,1) = linspace(startDate_m,endDate_m,36);
rainfallm_2012(:,2) = 1:1:12;

%% Daily rainfall plot
figure(1);clf;
subplot(1,2,1);hold on
boxplot(rain_avg)
plot(rainfallm_2012(:,2),rainfallm_2012(:,1),'-md','MarkerSize',8,'linewidth',1.5);
plot(rainfallm_2012(:,2),rainfallm_2013(:,1),'-gs','MarkerSize',8,'linewidth',1.5);
plot(rainfallm_2012(:,2),rainfallm_2014(:,1),'-bo','MarkerSize',8,'linewidth',1.5);
% axis equal square
axis tight;
xlabel('Months')
ylabel('Rainfall (mm)')
title('Monthly rainfall');
legend('2012','2013','2014')
months = {'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';'Jan';'Feb';'Mar'};
set(gca, 'XTick', 1:12); % center x-axis ticks on bins
set(gca, 'XTickLabel', months); % set x-axis labels
% grid on; grid minor,':';

subplot(1,2,2)
bar((1:1:36)',P_ETm);
axis tight
box on
xlabel('Days')
ylabel('Water Budget (mm)')
title('Weekly P-ET');
% xtick = {' ';' ';'Jun12';' ';' ';'Sep12';' ';' ';'Dec12';' ';' ';'Mar13';' ';' ';'Jun13';' ';' ';'Sep13';' ';' ';...
%     'Dec13';' ';' ';'Mar14';' ';' ';'Jun14';' ';' ';'Sep14';' ';' ';'Dec14';' ';' ';'Mar15'};
% set(gca, 'XTick', 1:36); % center x-axis ticks on bins
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', 3:3:36); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels
% datetick('x',12,'keepticks','keeplimits')

%% Drip data analysis
[D1 text1] = xlsread('2012-15_Stalagmate_data_golgotha.xlsx'); % Import Site 1 Drip data
text1(2,:) = [];
time1 = cell(size(text1,1),size(text1,2)/3);
stal_ID_1 = cell(1,size(text1,2)/3);
for i = 1:size(text1,2)/3
    time1(:,i) = text1(:,1+(i-1)*3);
    stal_ID_1(1,i) = text1(1,1+(i-1)*3);
end
time1 = time1(2:end,1:14);

[Drip2 text2] = xlsread('2012-15 Site 2_SWWA drip logger data.xlsx'); % Import Site 2 Drip data
% [Drip text] = xlsread('2012-15_Stalagmate_data_golgotha.xlsx'); % Import Site 1 Drip data
text2(2,:) = [];
time2 = cell(size(text2,1),size(text2,2)/3);
stal_ID_2 = cell(1,size(text2,2)/3);
for i = 1:size(text2,2)/3
    time2(:,i) = text2(:,1+(i-1)*3);
    stal_ID_2(1,i) = text2(1,1+(i-1)*3);
end
% time = time(1:82217,1:14); % Take the drip data upto 31/12/2014 as rainfall Site 1 data
time2 = time2(1:102286,:); % Take the drip data upto 31/12/2014 as rainfall Site 2 data

no_of_stalagmates2 = size(time2,2);
fprintf('Total Number of Stalagmates = %0.0f \n', no_of_stalagmates2)

i = 3:3:3*no_of_stalagmates2-1;
Drip2(:,i) = [];
Drip2 = Drip2(1:102286,:);

Data2 = zeros(size(Drip2,1),1);
Data2(1,1) = 0; % Site 1
for i = 1:size(Data2,1)-1
    Data2(i+1,1) = Data2(i,1)+(0.25/24);
end
Data2(:) = Data2(:) + datenum('05-Apr-2012 00:14:00');

for j = 1:size(Drip2,2)/2
    ind = find(isnan(Drip2(:,2*j)));
    no_data_location2{1,j} = 0;
    no_data_location2{1,j}(1,1) = Data2(ind(1),1);
    k = 0;
    for i = 1:size(ind,1)-1
        if ind(i+1) ~= ind(i)+1;
            k = k+2;
            no_data_location2{1,j}(k,1) = Data2(ind(i),1);
            no_data_location2{1,j}(k+1,1) = Data2(ind(i+1),1);
        end
    end
    no_data_location2{1,j}(k+2,1) = Data2(ind(i),1);
end

datestr(Data2(31872),'dd-mmm-yyyy HH:MM:SS');
datestr(Data2(66913),'dd-mmm-yyyy HH:MM:SS');
datestr(Data2(end),'dd-mmm-yyyy HH:MM:SS');
Data2_2012 = Data2(1:31872,:);
Data2_2013 = Data2(31873:66913,:);
Data2_2014 = Data2(66913:end,:);

%%
% [Drip_manual text_manual] = xlsread('Golgotha_manual_drips.per.min_kashif.xls'); % Import Manual Drip data
[Drip_manual text_manual] = xlsread('Golgotha_manual_drips.per.min.xls'); % Import Manual Drip data
formatIn = 'dd/mm/yyyy';
Data_manual_tick = datenum(text_manual(2:end,1),formatIn);
Data_manual = Data_manual_tick;
% Data_manual = Data_manual_tick(:) - Data_manual_tick(1) + 24;
Drip_manual = Drip_manual(:,1:5)*15; % Convert to 15 minutes drip data
% Drip_manual_norm(:,1) = mat2gray(Drip_manual(:,1));
% Drip_manual_norm(:,1) = Drip_manual(:,1)/min(Drip_manual(:,1));
% Drip_manual_norm(:,2) = Drip_manual(:,2)/min(Drip_manual(:,2));
% Drip_manual_norm(:,3) = Drip_manual(:,3)/min(Drip_manual(:,3));
% Drip_manual_norm(:,4) = Drip_manual(:,4)/min(Drip_manual(:,4));
% Drip_manual_norm(:,5) = Drip_manual(:,1)/min(Drip_manual(:,5));

[data_cumsum text_cumsum] = xlsread('cumsumprecip_FEW.xlsx'); % Import (CumSum Precipitation - FEW) data
time_cumsum = datenum(text_cumsum(249:355,1),formatIn);
cump_few = data_cumsum(249:355,1);
% cump_few(:) = cump_few(:)/min(cump_few);
cump_few_norm = (cump_few - min(cump_few)) / ( max(cump_few) - min(cump_few) );
Drip_manual_norm = zeros(size(Drip_manual));
for i = 1:size(Drip_manual,2)
    Drip_manual_norm(:,i) = (Drip_manual(:,i) - min(Drip_manual(:,i))) / ( max(Drip_manual(:,i)) - min(Drip_manual(:,i)) );
end 

%%
cump_few_movavg = MovAvg(cump_few_norm,10);
Drip_manual_movavg = zeros(size(Drip_manual));
for i = 1:size(Drip_manual,2)
    Drip_manual_movavg(:,i) = MovAvg(Drip_manual_norm(:,i),10);
end

%% Trend analysis of manual drip data and recharge/Cumulative P-ET anomaly
figure(2);clf;hold on
% subplot(1,2,1);hold on
% plot(Data_manual,Drip_manual_norm(:,1),'ro'); % Site 1A Manual Data
% plot(Data_manual,Drip_manual_norm(:,2),'go'); % Site 1B Manual Data
% plot(Data_manual,Drip_manual_norm(:,3),'bx'); % Site 2A Manual Data
% plot(Data_manual,Drip_manual_norm(:,4),'yx'); % Site 2B Manual Data
% plot(Data_manual,Drip_manual_norm(:,5),'cx'); % Site 2E Manual Data
% plot(time_cumsum,cump_few_norm,'md','MarkerSize',8);
% plot(Recharge_time,Recharge_norm,'k.');

% plot(Data_manual,Drip_manual_movavg(:,1),'ro'); % Site 1A Manual Data
% plot(Data_manual(19:end),Drip_manual_movavg(19:end,2),'go'); % Site 1B Manual Data
% plot(Data_manual(18:end),Drip_manual_movavg(18:end,3),'bx'); % Site 2A Manual Data
% plot(Data_manual,Drip_manual_movavg(:,4),'yx'); % Site 2B Manual Data
% plot(Data_manual(27:end),Drip_manual_movavg(27:end,5),'cx'); % Site 2E Manual Data
% plot(time_cumsum,cump_few_movavg,'md','MarkerSize',8);
% plot(Recharge_time,Recharge_norm,'k.');
% legend('1A Manual data','1B Manual data','2A Manual data','2B Manual data',...
%     '2E Manual data','Monthly CumSum(Precip - FEW) variation','Annual Recharge','Location','NorthWest')

% Drip_manual_smoothed = zeros(100*size(Data_manual,1),size(Drip_manual,2));
Data_manual_smoothed = linspace(Data_manual(1,1), Data_manual(end,1), 5000);
smoothedY = spline(Data_manual, Drip_manual_movavg(:,1), Data_manual_smoothed);
plot(Data_manual_smoothed,smoothedY,'r');

Data_manual_smoothed = linspace(Data_manual(19,1), Data_manual(end,1), 4000);
smoothedY = spline(Data_manual(19:end), Drip_manual_movavg(19:end,2), Data_manual_smoothed);
plot(Data_manual_smoothed,smoothedY,'g');

Data_manual_smoothed = linspace(Data_manual(18,1), Data_manual(end,1), 5000);
smoothedY = spline(Data_manual(18:end), Drip_manual_movavg(18:end,3), Data_manual_smoothed);
plot(Data_manual_smoothed,smoothedY,'b');

Data_manual_smoothed = linspace(Data_manual(2,1), Data_manual(end,1), 5000);
smoothedY = spline(Data_manual(2:end), Drip_manual_movavg(2:end,4), Data_manual_smoothed);
plot(Data_manual_smoothed,smoothedY,'y');

Data_manual_smoothed = linspace(Data_manual(27,1), Data_manual(end-1,1), 3500);
smoothedY = spline(Data_manual(27:end-1), Drip_manual_movavg(27:end-1,5), Data_manual_smoothed);
plot(Data_manual_smoothed,smoothedY,'c');

cumsum_smoothed = linspace(time_cumsum(1,1), time_cumsum(end,1), 5000);
smoothedY = spline(time_cumsum,cump_few_movavg, cumsum_smoothed);
plot(cumsum_smoothed,smoothedY,'m','LineWidth',3);

Recharge_smoothed = linspace(Recharge_time(1,1), Recharge_time(end,1), 1000);
smoothedY = spline(Recharge_time, Recharge_norm, Recharge_smoothed);
plot(Recharge_smoothed,smoothedY,'k','LineWidth',3);

legend('1A Manual data','1B Manual data','2A Manual data','2B Manual data',...
    '2E Manual data','Mean monthly anomaly of cumulative (P-ET) sum','Annual Recharge','Location','NorthWest')
% plot(period,DD,'--k*');
% axis equal square
axis tight
% xlim([Data(1,1) Data(end,1)])
xlabel('Days')
ylabel('Drips per 15 minutes (Normalized [0 1])') 
title('Comparison between Manual Drip rates and Variations of CUMSUM(Precip - ET) from monthly mean');
datetick('x',20,'keepticks','keeplimits')

%% Select Better quality Drip data
% Drip2(:,[11:12 19:20 27:28]) = []; % Choose only Better quality Drip data
% time2(:,[6 10 14]) = []; % Choose only Better quality Drip data time
Drip2(:,[9:10 29:30]) = []; % Choose only Better quality Drip data
time2(:,[5 15]) = []; % Choose only Better quality Drip data time
stal_ID_2(:,[5 15]) = []; % Choose only Better quality Drip data IDs
loc = Drip2(:,20)<= 2.5; % Discard a portion from site 2vii Drip data
Drip2(loc,19:20) = NaN;
% fprintf(' %18.58f\n',(Data2(35274)-Data2(35273)));
Drip2(find(Data2(:,1)==no_data_location2{1,1}(3,1)) : find(abs(Data2(:,1)-735480)<0.005),1:2) = NaN; % Discard a portion from site 2A Drip data
Drip2(find(abs(Data2(:,1)-735300)<0.005) : find(abs(Data2(:,1)-735350)<0.005),1:2) = NaN; % Discard a portion from site 2A Drip data


%% Plot 3 different hydrological years (2012, 2013 and 2014)
time_2012 = time2(1:31872,:);
time_2013 = time2(31873:66913,:);
time_2014 = time2(66913:end,:);
Drip_2012 = Drip2(1:31872,:);
Drip_2013 = Drip2(31873:66913,:);
Drip_2014 = Drip2(66913:end,:);

%% Import site 1 raw field (unfilled) data
[Data1_unfilled] = xlsread('Data1_unfilled.xlsx');
[Drip1_unfilled] = xlsread('Drip1_unfilled.xlsx');
stal_ID_1(:,[6 10 14]) = []; % Choose only Better quality Drip data IDs
[Data1] = xlsread('Data1.xlsx');
[Drip1] = xlsread('Drip1.xlsx');

%% Import Elevation data
[elv elv_text] = xlsread('Depth_data_site1_2.xlsx'); % Import elevation data
elv11 = elv(1:11,3);
elv1(1:2,1) = elv11(2:3,1); elv1(3,1) = elv11(1,1); elv1(4:11,1) = elv11(4:11,1);
elv22 = elv(1:18,9);
elv2(1:5,1) = elv22(2:6,1); elv2(6,1) = elv22(1,1); elv2(7:18,1) = elv22(7:18,1);

elevation_data(1:size(elv1,1),1) = elv1;
elevation_data(1+size(elv1,1):size(elv1,1)+size(elv2),1) = elv2;
depth11 = elv(1:11,4);
depth1(1:2,1) = depth11(2:3,1); depth1(3,1) = depth11(1,1); depth1(4:11,1) = depth11(4:11,1);
depth22 = elv(1:18,10);
depth2(1:5,1) = depth22(2:6,1); depth2(6,1) = depth22(1,1); depth2(7:18,1) = depth22(7:18,1);

depth_data(1:size(depth1,1),1) = depth1;
depth_data(1+size(depth1,1):size(depth1,1)+size(depth2),1) = depth2;

%% Drip data stat analysis
% Site 1
no_of_stalagmates1 = size(Drip1_unfilled,2)/2;
fprintf('Total Number of Stalagmates = %0.0f \n', no_of_stalagmates1)
avg_drip_1 = zeros(no_of_stalagmates1,1);
Variance_1 = zeros(no_of_stalagmates1,1);
std_1 = zeros(no_of_stalagmates1,1);
cov_1 = zeros(no_of_stalagmates1,1);
sk_1 = zeros(no_of_stalagmates1,1);
M_1 = zeros(no_of_stalagmates1,1);
figure(11);clf;figure(12);clf;
for j=1:no_of_stalagmates1
    Data1_unfilled(:,2) = Drip1_unfilled(:,1+(j-1)*2);
    Data1_unfilled(:,3) = Drip1_unfilled(:,2+(j-1)*2);
    Data_modified_1 = Data1_unfilled;
    % total_drip = sum(Data(:,3));
    Data_modified_1(isnan(Data_modified_1(:,2)),:) = [];
    
    total_drip_1 = sum(Data_modified_1(:,2));
    total_time_1 = size(Data_modified_1,1)/4;
    % total_time = Data_modified(end,1);
    avg_drip_1(j) = total_drip_1/total_time_1/4;
%     fprintf('Average drip rate per 15 minutes for %s Stalagmate = %0.1f \n', stal_ID_1{1,j}, avg_drip_1(j))
    
    Variance_1(j) = var(Data_modified_1(:,2));
%     fprintf('Variance in 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_1{1,j}, Variance_1(j))
    std_1(j) = std(Data_modified_1(:,2));
    cov_1(j) = std_1(j)/avg_drip_1(j)*100;
    
    figure(11);hold on
    subplot(4,5,j)
    hist(Data_modified_1(:,2),100)
    title(['Stalagmate ', stal_ID_1{1,j}])
%     title(['TI = ', num2str(t) ' with Variable = ', num2str(j)]);
    sk_1(j) = skewness(Data_modified_1(:,2));
%     fprintf('Skewness of 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_1{1,j}, sk_1(j))
    M_1(j) = mode(Data_modified_1(:,2));
%     fprintf('Mode of 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_1{1,j}, M_1(j))
    
%     figure(12);hold on
%     subplot(4,5,j)
%     autocorr(Data_modified_1(:,2),1000)
%     title(['Stalagmate ', stal_ID_1{1,j}])
end

% Site 2
no_of_stalagmates2 = size(time2,2);
fprintf('Total Number of Stalagmates = %0.0f \n', no_of_stalagmates2)
avg_drip_2 = zeros(no_of_stalagmates2,1);
Variance_2 = zeros(no_of_stalagmates2,1);
std_2 = zeros(no_of_stalagmates2,1);
cov_2 = zeros(no_of_stalagmates2,1);
sk_2 = zeros(no_of_stalagmates2,1);
M_2 = zeros(no_of_stalagmates2,1);
figure(13);clf;figure(14);clf;
for j=1:no_of_stalagmates2
    Data2(:,2) = Drip2(:,1+(j-1)*2);
    Data2(:,3) = Drip2(:,2+(j-1)*2);
    Data_modified_2 = Data2;
    % total_drip = sum(Data(:,3));
    Data_modified_2(isnan(Data_modified_2(:,2)),:) = [];
    
    total_drip_2 = sum(Data_modified_2(:,2));
    total_time_2 = size(Data_modified_2,1)/4;
    % total_time = Data_modified(end,1);
    avg_drip_2(j) = total_drip_2/total_time_2/4;
%     fprintf('Average drip rate per 15 minutes for %s Stalagmate = %0.1f \n', stal_ID_2{1,j}, avg_drip_2(j))
    
    Variance_2(j) = var(Data_modified_2(:,2));
%     fprintf('Variance in 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_2{1,j}, Variance_2(j))
    std_2(j) = std(Data_modified_2(:,2));
    cov_2(j) = std_2(j)/avg_drip_2(j)*100;
    
    figure(13);hold on
    subplot(4,5,j)
    hist(Data_modified_2(:,2),100)
    title(['Stalagmate ', stal_ID_2{1,j}])
%     title(['TI = ', num2str(t) ' with Variable = ', num2str(j)]);
    sk_2(j) = skewness(Data_modified_2(:,2));
%     fprintf('Skewness of 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_2{1,j}, sk_2(j))
    M_2(j) = mode(Data_modified_2(:,2));
%     fprintf('Mode of 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_2{1,j}, M_2(j))
    
%     figure(14);hold on
%     subplot(4,5,j)
%     autocorr(Data_modified_2(:,2),1000)
%     title(['Stalagmate ', stal_ID_2{1,j}])
end

avg_drip(1:size(avg_drip_1,1),1) = avg_drip_1;
avg_drip(1+size(avg_drip_1,1):size(avg_drip_1,1)+size(avg_drip_2,1),1) = avg_drip_2;
cov(1:size(cov_1,1),1) = cov_1;
cov(1+size(cov_1,1):size(cov_1,1)+size(cov_2,1),1) = cov_2;
sk(1:size(sk_1,1),1) = sk_1;
sk(1+size(sk_1,1):size(sk_1,1)+size(sk_2,1),1) = sk_2;
M(1:size(M_1,1),1) = M_1;
M(1+size(M_1,1):size(M_1,1)+size(M_2,1),1) = M_2;

% flow_class = [1 1 1 1 1 1 4 1 3 2 1 1 1 3 2 4 1 2 3 2 1 3 1 4 1 1 1 2 1];
flow_class = [2 2 2 2 2 2 1 2 3 4 2 2 2 3 4 1 2 4 3 4 2 3 2 1 2 2 2 4 2];

%% Limestone thickness plots
figure(3);clf;hold on;
% subplot(2,2,1);
% scatter(elevation_data,avg_drip,1,'filled');
scatter(log(avg_drip_1),depth1,20,flow_class(1:size(avg_drip_1,1))','d','LineWidth',1.5);
scatter(log(avg_drip_2),depth2,20,flow_class(size(avg_drip_1,1)+1:end)','o','LineWidth',1.5);
% scatter(avg_drip,depth_data,50,flow_class,'LineWidth',2);
% scatter(avg_drip_1,depth1,20,6*ones(size(avg_drip_1,1),1),'filled','d');
% scatter(avg_drip_2,depth2,20,5*ones(size(avg_drip_2,1),1),'filled','o');
% scatter(avg_drip,depth_data,50,flow_class,'LineWidth',2);
axis equal square
axis tight; box on
% xlim([0 60])
ylabel('Limestone thickness (m)')
xlabel('log (Avg Drip per 15 mins)')
legend('Site 1','Site 2')
% legend('Site 1','Site 2','Flow type')
title('Drip rate vs Limestone thickness');
% r2_1 = (corr(log(avg_drip_1),elv1))^2;
% r2_2 = (corr(log(avg_drip_2),elv2))^2;
r2 = (corr(avg_drip,elevation_data))^2;

figure(4);clf;hold on;
% subplot(2,2,1);
% scatter(elevation_data,avg_drip,1,'filled');
scatter(sk_1,depth1,20,flow_class(1:size(avg_drip_1,1))','d','LineWidth',1.5);
scatter(sk_2,depth2,20,flow_class(size(avg_drip_1,1)+1:end)','o','LineWidth',1.5);
% scatter(sk_1,depth1,20,6*ones(size(sk_1,1),1),'filled','d');
% scatter(sk_2,depth2,20,5*ones(size(sk_2,1),1),'filled','o');
% scatter(sk,depth_data,50,flow_class,'LineWidth',2);
axis equal square
axis tight; box on
% xlim([0 60])
ylabel('Limestone thickness (m)')
xlabel('Skewness')
legend('Site 1','Site 2')
title('Skewness vs Limestone thickness');
r3 = (corr(sk,elevation_data))^2;

figure(5);clf;hold on;
% subplot(2,2,1);
% scatter(elevation_data,avg_drip,1,'filled');
scatter(log(cov_1),depth1,20,flow_class(1:size(avg_drip_1,1))','d','LineWidth',1.5);
scatter(log(cov_2),depth2,20,flow_class(size(avg_drip_1,1)+1:end)','o','LineWidth',1.5);
% scatter(cov_1,depth1,20,6*ones(size(cov_1,1),1),'filled','d');
% scatter(cov_2,depth2,20,5*ones(size(cov_2,1),1),'filled','o');
% scatter(cov,depth_data,50,flow_class,'LineWidth',2);
axis equal square
axis tight; box on
% xlim([0 60])
ylabel('Limestone thickness (m)')
xlabel('log (COV)')
legend('Site 1','Site 2')
title('Coefficient of Variation (COV) vs Limestone thickness');
r4 = (corr(cov,elevation_data))^2;


%% Find locations of missing data
for j = 1:size(Drip2,2)/2
    ind = find(isnan(Drip2(:,2*j)));
    no_data_location2{1,j} = 0;
    no_data_location2{1,j}(1,1) = Data2(ind(1),1);
    k = 0;
    for i = 1:size(ind,1)-1
        if ind(i+1) ~= ind(i)+1;
            k = k+2;
            no_data_location2{1,j}(k,1) = Data2(ind(i),1);
            no_data_location2{1,j}(k+1,1) = Data2(ind(i+1),1);
        end
    end
    no_data_location2{1,j}(k+2,1) = Data2(ind(i+1),1);
end

%% Fill data gaps with interpolated values
for j = 1:size(Drip2,2)/2
    x = Drip2(:,2*j);
    movavg = MovAvg(x,500);
%     figure(50+j);clf;hold on;
%     plot(Data2(:,1),x,':y'); % Site 1A Data quality = Good
%     plot(Data2(:,1),movavg,':b'); % Site 1i Data quality = Good
%     axis equal square
%     axis tight
%     title(['Drip data time series for Site ', stal_ID_2{1,j}])
%     xlabel('Days')
%     ylabel('Drips per 15 minutes')
    
    cord = zeros(size(no_data_location2{1,j},1),1);
    for i = 1:size(no_data_location2{1,j},1)
        if bitget(i,1)
            cord(i) = no_data_location2{1,j}(i) - datenum('00-000-0000 12:00:00');
        else
            cord(i) = no_data_location2{1,j}(i) + datenum('00-000-0000 12:00:00');
            movavg ( find(Data2==cord(i-1)) : find(Data2==cord(i)) ) = nan;
        end
    end
    % ind = find(isnan(movavg));
    
    nans = isnan(movavg);
    movavg(nans) = interp1(Data2(~nans), movavg(~nans), Data2(nans), 'pchip'); % linear spline pchip
    plot(Data2,movavg,'--c');
    if no_data_location2{1,j}(1) == Data2(1)
        no_data_location2{1,j}(1:2) = [];
    end
    for i = 0:size(no_data_location2{1,j},1)/2-1
        plot(no_data_location2{1,j}(1+2*i,1),Drip2(find(Data2(:,1)==no_data_location2{1,j}(1+2*i,1))-1,2*j),'rx','LineWidth',3,'MarkerSize',10);
    end
    legend('Drip data','Mov Avg','Interpolated data','Location of missing data','Location','Best')
    datetick('x',20,'keepticks','keeplimits')
    
%     y = Drip(:,2*j);
    if no_data_location2{1,j}(end) == Data2(end)
        no_data_location2{1,j}(end-1:end) = [];
    end
    for i = 0:size(no_data_location2{1,j},1)/2-1  % Fill the gaps in original drip data
        x ( find(Data2==no_data_location2{1,j}(2*i+1)-datenum('00-000-0000 12:00:00')) : find(Data2==no_data_location2{1,j}(2*i+2)+datenum('00-000-0000 12:00:00')) ) = movavg ( find(Data2==no_data_location2{1,j}(2*i+1)-datenum('00-000-0000 12:00:00')) : find(Data2==no_data_location2{1,j}(2*i+2)+datenum('00-000-0000 12:00:00')) );
        Drip2(:,2*j) = x;
    end
end

%% Plot soda-straw drip data
figure(6);clf;hold on
plot(Data1(:,1),Drip1(:,14),'--k'); % Site 1v Data quality = Good
plot(Data2(:,1),Drip2(:,10),'--r'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,26),'--g'); % Site 2xi Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('1v','2iii','2xi')
title('Soda-straw drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

% Plot site 1 icicle flow drip data
figure(7);clf;hold on
plot(Data1(:,1),Drip1(:,2),'--k'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,4),'--r'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,6),'--g'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,8),'--b'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,10),'--m'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,12),'--y'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,16),':k'); % Site 1v Data quality = Good
plot(Data1(:,1),Drip1(:,22),':r'); % Site 1v Data quality = Good
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('1A','1B','1i','1ii','1iii','1ix','1vi','1xi')
title('Site 1 icicle flow drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

% Plot site 2 icicle flow drip data
figure(8);clf;hold on
plot(Data2(:,1),Drip2(:,2),'--k'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,4),'--r'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,12),'--g'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,20),'--b'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,24),'--m'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,28),'--y'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,30),':k'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,32),':g'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,36),':r'); % Site 2iii Data quality = Good
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('2A','2B','2vi','2vii','2x','2xii','2xiv','2xv','2xvii')
title('Site 2 icicle flow drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

% Plot combined flow drip data
figure(9);clf;hold on
plot(Data1(:,1),Drip1(:,18),'--k'); % Site 1v Data quality = Good
plot(Data2(:,1),Drip2(:,6),'--r'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,16),'--g'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,22),'--g'); % Site 2xi Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('1viii','2E','2v','2viii')
title('Combined flow drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

% Plot fracture flow drip data
figure(10);clf;hold on
plot(Data1(:,1),Drip1(:,20),'--k'); % Site 1v Data quality = Good
plot(Data2(:,1),Drip2(:,8),'--r'); % Site 2iii Data quality = Good
plot(Data2(:,1),Drip2(:,14),'--g'); % Site 2xi Data quality = Better
% plot(Data2(:,1),Drip2(:,18),'--b'); % Site 2xi Data quality = Better
plot(Data2(:,1),Drip2(:,34),'--m'); % Site 2xi Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('1x','2i','2ix','2xvi')
title('Fracture flow drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels


% Plot fracture flow drip data
figure(11);clf;hold on
plot(Data2(:,1),Drip2(:,18),'--b'); % Site 2xi Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('2vi')
title('Fracture flow drip sites');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

%% Plot filled final drip data
figure(15);clf;
subplot(1,2,1);hold on
plot(Data2(:,1),Drip2(:,2),'--k'); % Site 2B Data quality = Good
plot(Data2(:,1),Drip2(:,4),'--m'); % Site 2iii Data quality = Better
plot(Data2(:,1),Drip2(:,10),'--r'); % Site 2iv Data quality = Good
% plot(Data2(:,1),Drip2(:,12),'--g'); % Site 2iv Data quality = Good
% plot(Data2(:,1),Drip2(:,20),'--b'); % Site 2xi Data quality = Best
% plot(Data2(:,1),Drip2(:,24),':r'); % Site 2xiv Data quality = Better
plot(Data2(:,1),Drip2(:,26),':g'); % Site 2xv Data quality = Better
plot(Data2(:,1),Drip2(:,28),':c'); % Site 2xvii Data quality = Better
plot(Data2(:,1),Drip2(:,30),':b'); % Site 2xiv Data quality = Better
% plot(Data2(:,1),Drip2(:,32),':k'); % Site 2xv Data quality = Better
% plot(Data2(:,1),Drip2(:,36),':m'); % Site 2xvii Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('2A','2B','2iii','2xi','2xiii','2xiv')
% legend('2A','2B','2iii','2iv','2vii','2x','2xi','2xiii','2xiv','2xv','2xvii')
title('Daily drip rates from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

subplot(1,2,2);hold on;
bar(rainfall(5:end,2),P_ETw);
plot(Recharge_smoothed,smoothedY*100,'m','LineWidth',3);
axis equal square
axis tight
box on
xlim([Data2(1,1) Data2(end,1)])
ylim([0 120])
xlabel('Days')
ylabel('Rainfall (mm)')
legend('Weekly (P-ET)','Estimated recharge trend')
title('Weekly P-ET from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(rainfall(12,2),rainfall(end,2),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

figure(16);clf;
subplot(1,2,1);hold on
plot(Data2(:,1),Drip2(:,6),'--k'); % Site 2E Data quality = Good
plot(Data2(:,1),Drip2(:,8),'--g'); % Site 2iv Data quality = Good
plot(Data2(:,1),Drip2(:,14),'--m'); % Site 2i Data quality = Better
plot(Data2(:,1),Drip2(:,16),'--y'); % Site 2x Data quality = Best
plot(Data2(:,1),Drip2(:,22),'--r'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2(:,34),'--c'); % Site 2xvii Data quality = Better
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('2E','2i','2ix','2v','2viii','2xvi')
title('Daily drip rates from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

subplot(1,2,2);hold on;
bar(rainfall(5:end,2),P_ETw);
plot(Recharge_smoothed,smoothedY*100,'m','LineWidth',3);
axis equal square
axis tight
box on
xlim([Data2(1,1) Data2(end,1)])
xlabel('Days')
ylabel('Rainfall (mm)')
title('Weekly P-ET from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(rainfall(12,2),rainfall(end,2),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

figure(17);clf;
subplot(1,2,1);hold on
plot(Data2(:,1),Drip2(:,18),':r'); % Site 2ix Data quality = Good
axis equal square
axis tight
xlabel('Days') 
ylabel('Drips per 15 minutes') 
legend('2vi')
title('Daily drip rates from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(Data2(7000,1),Data2(end,1),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

subplot(1,2,2);hold on;
bar(rainfall(5:end,2),P_ETw);
plot(Recharge_smoothed,smoothedY*100,'m','LineWidth',3);
axis equal square
axis tight
box on
xlim([Data2(1,1) Data2(end,1)])
ylim([min(P_ETw)-1 120])
xlabel('Days')
ylabel('Rainfall (mm)')
title('Weekly P-ET from 05-04-2012 to 06-03-2015');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun12';'Sep12';'Dec12';'Mar13';'Jun13';'Sep13';'Dec13';'Mar14';'Jun14';'Sep14';'Dec14';'Mar15'};
set(gca, 'XTick', linspace(rainfall(12,2),rainfall(end,2),12)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

%% Test drip response to certain rainfall events
% [Data1] = xlsread('Data1.xlsx');
% [Drip1] = xlsread('Drip1.xlsx');
Drip1_movavg(:,18) = MovAvg(Drip1(:,18),2000);
Drip1_movavg(:,20) = MovAvg(Drip1(:,20),2000);
Drip2_movavg(:,6) = MovAvg(Drip2(:,6),2000);
Drip2_movavg(:,8) = MovAvg(Drip2(:,8),2000);
Drip2_movavg(:,14) = MovAvg(Drip2(:,14),2000);
Drip2_movavg(:,16) = MovAvg(Drip2(:,16),2000);
Drip2_movavg(:,22) = MovAvg(Drip2(:,22),2000);
Drip2_movavg(:,34) = MovAvg(Drip2(:,34),2000);
Drip2_movavg(:,4) = MovAvg(Drip2(:,4),2000);
Drip2_movavg(:,28) = MovAvg(Drip2(:,28),2000);

%% Plot
figure(18);clf;
subplot(1,3,1);hold on
bar(rainfall(5:end,2),P_ETw/5,0.6);
% plot(Data1(:,1),Drip1(:,2),'--k'); % Site 2E Data quality = Good
% plot(Data1(:,1),Drip1(:,4),'--g'); % Site 2iv Data quality = Good
% plot(Data1(:,1),Drip1(:,6),'--m'); % Site 2i Data quality = Better
% plot(Data1(:,1),Drip1(:,8),'--y'); % Site 2x Data quality = Best
% plot(Data1(:,1),Drip1(:,10),'--k'); % Site 2E Data quality = Good
% plot(Data1(:,1),Drip1(:,12),'--g'); % Site 2iv Data quality = Good
% plot(Data1(:,1),Drip1(:,14),'--m'); % Site 2i Data quality = Better
% plot(Data1(:,1),Drip1(:,16),'--y'); % Site 2x Data quality = Best
plot(Data1(:,1),Drip1(:,18),':r'); % Site 2xiii Data quality = Better
% plot(Data1(:,1),Drip1(:,20),':g'); % Site 2xvii Data quality = Better
plot(Data1(:,1),Drip1_movavg(:,18),'c'); % Site 2xiii Data quality = Better
% plot(Data1(:,1),Drip1_movavg(:,20),'k'); % Site 2xvii Data quality = Better
legend('P-ET','1viii','1viii trend')
axis equal square;axis tight;box on
xlim([datenum('03/03/2013',formatIn) datenum('02/03/2014',formatIn)])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates for site 1');
datetick('x',20,'keepticks','keeplimits')

subplot(1,3,2);hold on;
bar(rainfall(5:end,2),P_ETw,0.6);
% plot(Data2(:,1),Drip2(:,6),':k'); % Site 2E Data quality = Good
% plot(Data2(:,1),Drip2(:,8),':g'); % Site 2iv Data quality = Good
plot(Data2(:,1),Drip2(:,14),':m'); % Site 2i Data quality = Better
% plot(Data2(:,1),Drip2(:,16),':y'); % Site 2x Data quality = Best
plot(Data2(:,1),Drip2(:,22),':g'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2(:,34),':r'); % Site 2xvii Data quality = Better
% legend('P-ET','2E','2i','2ix','2v','2viii','2xvi')
% plot(Data2(:,1),Drip2_movavg(:,6),'c'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2_movavg(:,8),'c'); % Site 2xvii Data quality = Better
plot(Data2(:,1),Drip2_movavg(:,14),'c'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2_movavg(:,16),'c'); % Site 2xvii Data quality = Better
plot(Data2(:,1),Drip2_movavg(:,22),'y'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2_movavg(:,34),'k'); % Site 2xvii Data quality = Better
legend('P-ET','2ix','2viii','2xvi','2ix trend','2viii trend','2xvi trend')
axis equal square;axis tight;box on
xlim([datenum('03/03/2013',formatIn) datenum('02/03/2014',formatIn)])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates (intermediate) for site 2');
datetick('x',20,'keepticks','keeplimits')


subplot(1,3,3);hold on;
bar(rainfall(5:end,2),P_ETw/10,0.6);
% plot(Data2(:,1),Drip2(:,2),':k'); % Site 2E Data quality = Good
plot(Data2(:,1),Drip2(:,4),':g'); % Site 2iv Data quality = Good
% plot(Data2(:,1),Drip2(:,12),':m'); % Site 2i Data quality = Better
% plot(Data2(:,1),Drip2(:,10),':y'); % Site 2x Data quality = Best
% plot(Data2(:,1),Drip2(:,20),':b'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2(:,24),':r'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,26),':c'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2(:,28),':r'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,30),'--g'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2(:,32),'--m'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,36),'--b'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2_movavg(:,4),'c'); % Site 2xiii Data quality = Better
plot(Data2(:,1),Drip2_movavg(:,28),'k'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2_movavg(:,14),'c'); % Site 2xiii Data quality = Better
% % plot(Data2(:,1),Drip2_movavg(:,16),'c'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2_movavg(:,22),'k'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2_movavg(:,34),'y'); % Site 2xvii Data quality = Better
legend('P-ET','2B','2xiii','2B trend','2xiii trend')
axis equal square;axis tight;box on
xlim([datenum('03/03/2013',formatIn) datenum('02/03/2014',formatIn)])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates (low) for site 2');
datetick('x',20,'keepticks','keeplimits')

%% Format Daily data and rainfall
formatIn = 'yyyy-mm-dd HH:MM:SS';
Data1_day = Data1 (find(abs(Data1-datenum('2013-03-03 00:00:00',formatIn))<0.005) : find(abs(Data1-datenum('2014-03-03 00:00:00',formatIn))<0.005),1);
Drip1_day = Drip1 (find(abs(Data1-datenum('2013-03-03 00:00:00',formatIn))<0.005) : find(abs(Data1-datenum('2014-03-03 00:00:00',formatIn))<0.005),:);
Drip2_day = Drip2 (find(abs(Data2-datenum('2013-03-03 00:00:00',formatIn))<0.005) : find(abs(Data2-datenum('2014-03-03 00:00:00',formatIn))<0.005),:);
Drip1_d = zeros(365,size(Drip1_day,2));
Drip2_d = zeros(365,size(Drip2_day,2));
Data1_d(1:365,1) = datenum('2013-03-03 00:00:00',formatIn):1:datenum('2014-03-02 00:00:00',formatIn);
Data2_d(1:365,1) = datenum('2013-03-03 00:00:00',formatIn):1:datenum('2014-03-02 00:00:00',formatIn);
for i = 1:365
    Drip1_d(i,:) = sum(Drip1_day(1+(i-1)*96:(i*96),:))/96;
    Drip2_d(i,:) = sum(Drip2_day(1+(i-1)*96:(i*96),:))/96;
end
rainfall_d = rainfalld(3349:3713,:);


%% Daily plots
figure(19);clf;
subplot(1,3,1);hold on
bar(Data1_d,rainfall_d/2.5,0.6);
% plot(Data1(:,1),Drip1(:,2),'--k'); % Site 2E Data quality = Good
% plot(Data1(:,1),Drip1(:,4),'--g'); % Site 2iv Data quality = Good
% plot(Data1(:,1),Drip1(:,6),'--m'); % Site 2i Data quality = Better
% plot(Data1(:,1),Drip1(:,8),'--y'); % Site 2x Data quality = Best
% plot(Data1(:,1),Drip1(:,10),'--k'); % Site 2E Data quality = Good
% plot(Data1(:,1),Drip1(:,12),'--g'); % Site 2iv Data quality = Good
% plot(Data1(:,1),Drip1(:,14),'--m'); % Site 2i Data quality = Better
% plot(Data1(:,1),Drip1(:,16),'--y'); % Site 2x Data quality = Best
plot(Data1_d(:,1),Drip1_d(:,18),'r'); % Site 2xiii Data quality = Better
% plot(Data1(:,1),Drip1(:,20),':g'); % Site 2xvii Data quality = Better
legend('Daily precipitation','1viii')
axis equal square;axis tight;box on
ylim([0 40])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates for site 1');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun13';'Sep13';'Dec13';'Mar14'};
set(gca, 'XTick', linspace(Data1_d(90,1),Data1_d(end,1),4)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

% figure(19);clf;hold on;
subplot(1,3,2);hold on;
bar(Data1_d,rainfall_d,0.6);
% plot(Data2(:,1),Drip2(:,6),':k'); % Site 2E Data quality = Good
% plot(Data2(:,1),Drip2(:,8),':g'); % Site 2iv Data quality = Good
% plot(Data2(:,1),Drip2(:,26)*100,':y'); % Site 2x Data quality = Best
plot(Data2_d(:,1),Drip2_d(:,14),'m'); % Site 2i Data quality = Better
plot(Data2_d(:,1),Drip2_d(:,22),'g'); % Site 2xiii Data quality = Better
plot(Data2_d(:,1),Drip2_d(:,34),'r'); % Site 2xvii Data quality = Better
% legend('P-ET','2E','2i','2ix','2v','2viii','2xvi')
legend('Daily precipitation','2ix','2viii','2xvi')
axis equal square;axis tight;box on
ylim([0 100])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates (intermediate) for site 2');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun13';'Sep13';'Dec13';'Mar14'};
set(gca, 'XTick', linspace(Data1_d(90,1),Data1_d(end,1),4)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

subplot(1,3,3);hold on;
bar(Data1_d,rainfall_d/5,0.6);
% plot(Data2(:,1),Drip2(:,2),':k'); % Site 2E Data quality = Good
plot(Data2_d(:,1),Drip2_d(:,4),'g'); % Site 2iv Data quality = Good
% plot(Data2(:,1),Drip2(:,12),':m'); % Site 2i Data quality = Better
% plot(Data2(:,1),Drip2(:,10),':y'); % Site 2x Data quality = Best
% plot(Data2(:,1),Drip2(:,20),':b'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2(:,24),':r'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,26),':c'); % Site 2xiii Data quality = Better
plot(Data2_d(:,1),Drip2_d(:,28),'r'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,30),'--g'); % Site 2xiii Data quality = Better
% plot(Data2(:,1),Drip2(:,32),'--m'); % Site 2xvii Data quality = Better
% plot(Data2(:,1),Drip2(:,36),'--b'); % Site 2xiii Data quality = Better
legend('Daily precipitation','2B','2xiii')
axis equal square;axis tight;box on
ylim([0 20])
xlabel('Days');ylabel('Drips per 15 minutes');title('Drip rates (low) for site 2');
% datetick('x',20,'keepticks','keeplimits')
xtick = {'Jun13';'Sep13';'Dec13';'Mar14'};
set(gca, 'XTick', linspace(Data1_d(90,1),Data1_d(end,1),4)); % center x-axis ticks on bins
set(gca, 'XTickLabel', xtick); % set x-axis labels

%% Plot 3 different hydrological years (2012, 2013 and 2014)
Drip_2012 = Drip2(1:31872,:);
Drip_2013 = Drip2(31873:66913,:);
Drip_2014 = Drip2(66913:end,:);


%% Validate automatic Stalagmates drip data with manual drip data
% [Data1] = xlsread('Data1.xlsx');
% [Drip1] = xlsread('Drip1.xlsx');

figure(20);clf;
subplot(2,2,1);hold on
plot(Data1(:,1),Drip1(:,2),':c'); % Site 1A Data quality = Good
plot(Data_manual,Drip_manual(:,1),'ks','LineWidth',1.5,'MarkerSize',4); % Site 1A Manual Data
plot(Data1(:,1),Drip1(:,4),':b'); % Site 1B Data quality = Better
plot(Data_manual,Drip_manual(:,2),'kd','LineWidth',1.5,'MarkerSize',5); % Site 1B Manual Data
% plot(Data2(:,1),Drip2(:,2),':r'); % Site 1A Data quality = Good
% plot(Data_manual,Drip_manual(:,3),'kx','LineWidth',1,'MarkerSize',5); % Site 2A Manual Data
% plot(Data2(:,1),Drip2(:,4),':m'); % Site 1B Data quality = Better
% plot(Data_manual,Drip_manual(:,4),'k+','LineWidth',1,'MarkerSize',5); % Site 2B Manual Data
% plot(Data2(:,1),Drip2(:,6),':g'); % Site 1A Data quality = Good
% plot(Data_manual,Drip_manual(:,5),'go','LineWidth',5,'MarkerSize',5); % Site 2E Manual Data
legend('1A Loggers data','1A Manual data','1B Loggers data','1B Manual data','Location','NorthWest')
axis tight
xlim([Data1(1,1) Data1(end,1)]); ylim([1.5 5.5]);
xlabel('Days')
ylabel('Drips per 15 minutes') 
title('Drip rates from 01-09-2012 to 06-03-2015');
datetick('x',20,'keepticks','keeplimits')

subplot(2,2,2);hold on
plot(Data2(:,1),Drip2(:,2),':r'); % Site 1A Data quality = Good
plot(Data_manual,Drip_manual(:,3),'kx','LineWidth',2,'MarkerSize',7); % Site 2A Manual Data
legend('2A Loggers data','2A Manual data','Location','NorthWest')
axis tight
xlim([Data1(1,1) Data1(end,1)]); ylim([0 4]);
xlabel('Days')
ylabel('Drips per 15 minutes') 
title('Drip rates from 01-09-2012 to 06-03-2015');
datetick('x',20,'keepticks','keeplimits')

subplot(2,2,3);hold on
plot(Data2(:,1),Drip2(:,4),':m'); % Site 1B Data quality = Better
plot(Data_manual,Drip_manual(:,4),'k+','LineWidth',2,'MarkerSize',6); % Site 2B Manual Data
legend('2B Loggers data','2B Manual data','Location','NorthWest')
axis tight
xlim([Data1(1,1) Data1(end,1)]); ylim([2 5]);
xlabel('Days')
ylabel('Drips per 15 minutes') 
title('Drip rates from 01-09-2012 to 06-03-2015');
datetick('x',20,'keepticks','keeplimits')

subplot(2,2,4);hold on
plot(Data2(:,1),Drip2(:,6),':g'); % Site 1A Data quality = Good
plot(Data_manual,Drip_manual(:,5),'k*','LineWidth',1,'MarkerSize',5); % Site 2E Manual Data
legend('2E Loggers data','2E Manual data','Location','NorthWest')
axis tight
xlim([Data1(1,1) Data1(end,1)]); ylim([21 32]);
xlabel('Days')
ylabel('Drips per 15 minutes') 
title('Drip rates from 01-09-2012 to 06-03-2015');
datetick('x',20,'keepticks','keeplimits')

%% Drip counts and water volume measurements for Stalagmates
drip_count_1 = cell(1,size(Drip1,2)/2);
for j=1:size(Drip1,2)/2
    drip_count_1{1,j} = Drip1(:,2*j-1);
end
for j=1:size(Drip1,2)/2
    drip_count_1{1,j}(isnan(drip_count_1{1,j}(:,1)),:) = [];
end
total_drip_count_1 = zeros(1,size(Drip1,2)/2);
for j=1:size(Drip1,2)/2
    total_drip_count_1(1,j) = sum(drip_count_1{1,j});
end
total_water_volume_1 = total_drip_count_1(1,:)*0.1433/1000;
t_1 = zeros(1,size(Drip1,2)/2);
for j=1:size(Drip1,2)/2
    t_1(1,j) = (size(drip_count_1{1,j},1)-1)*0.25/24/365;
end
%
drip_count_2 = cell(1,size(Drip2,2)/2);
for j=1:size(Drip2,2)/2
    drip_count_2{1,j} = Drip2(:,2*j-1);
end
for j=1:size(Drip2,2)/2
    drip_count_2{1,j}(isnan(drip_count_2{1,j}(:,1)),:) = [];
end
total_drip_count_2 = zeros(1,size(Drip2,2)/2);
for j=1:size(Drip2,2)/2
    total_drip_count_2(1,j) = sum(drip_count_2{1,j});
end
total_water_volume_2 = total_drip_count_2(1,:)*0.1433/1000;
t_2 = zeros(1,size(Drip2,2)/2);
for j=1:size(Drip2,2)/2
    t_2(1,j) = (size(drip_count_2{1,j},1)-1)*0.25/24/365;
end

%%
water_volume_1 = total_water_volume_1(1,:)./t_1(1,:);
water_volume_2 = total_water_volume_2(1,:)./t_2(1,:);
water_volume(1:size(water_volume_1,2),1) = water_volume_1;
water_volume(1+size(water_volume_1,2):size(water_volume_1,2)+size(water_volume_2,2),1) = water_volume_2;

for j=1:size(water_volume_1,2)
    fprintf('Total recharge volume for %s Stalagmate = %0.2f (L/yr) \n', stal_ID_1{1,j}, water_volume_1(1,j))
end
for j=1:size(water_volume_2,2)
    fprintf('Total recharge volume for %s Stalagmate = %0.2f (L/yr) \n', stal_ID_2{1,j}, water_volume_2(1,j))
end

total_rainfall = sum(rainfall(:,1));
avg_rainfall = total_rainfall/size(rainfall,1);
fprintf('Average rainfall = %0.2f (mm/week) \n', avg_rainfall)
rainfall_volume = avg_rainfall/1000*52*(1*1)*1000; % Unit = L/yr in 1 m2 area
fprintf('Total rainfall volume = %0.2f (L/yr/m2) \n', rainfall_volume)

% std_dis(j) = std(water_volume);
% cov_dis(j) = std_dis(j)/avg_drip_2(j)*100;
    
%% Sampling interval VS COV test
avg_drip_1 = zeros(no_of_stalagmates1,5);
% Variance_2 = zeros(no_of_stalagmates2,1);
std_1 = zeros(no_of_stalagmates1,5);
cov_1 = zeros(no_of_stalagmates1,5);
k = [1 4 96 672 2880];
figure(12);clf;figure(14);clf;figure(112);clf;figure(114);clf;
for i = 1:5
    for j = 1:no_of_stalagmates1
        Data1_unfilled(:,2) = Drip1_unfilled(:,1+(j-1)*2);
        %         Data1_unfilled(:,3) = Drip1_unfilled(:,2+(j-1)*2);
        Data_modified_1 = zeros(size(Data1_unfilled,1),2);
        Data_modified_1(:,1) = Data1_unfilled(:,1);
        Data_modified_1(:,2) = filter( ones(k(i),1), 1, Data1_unfilled(:,2) );
        Data_modified_1 = Data_modified_1(k(i):k(i):end,:);
        
%         Data_modified_1 = Data1_unfilled(1:k(i):end,:);
        
        Data_modified_1(isnan(Data_modified_1(:,2)),:) = [];
        
        total_drip_1 = sum(Data_modified_1(:,2));
        total_time_1 = size(Data_modified_1,1)/4*k(i); % Time in hours % DateNumber = datenum(total_time_1); DateString = datestr(total_time_1);
        % total_time = Data_modified(end,1);
        avg_drip_1(j,i) = total_drip_1/total_time_1/4;
        %     fprintf('Average drip rate per 15 minutes for %s Stalagmate = %0.1f \n', stal_ID_1{1,j}, avg_drip_1(j))
        
        %     Variance_1(j) = var(Data_modified_1(:,2));
        %     fprintf('Variance in 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_1{1,j}, Variance_1(j))
        std_1(j,i) = std(Data_modified_1(:,2));
%         cov_1(j,i) = std_1(j,i)/avg_drip_1(j,i)*100;
        cov_1(j,i) = std_1(j,i)/mean(Data_modified_1(:,2))*100;
        %ACF plots with daily frequency
        if i==3
            figure(12);hold on
            subplot(4,5,j)
            autocorr(Data_modified_1(:,2),200)
            title(['Stalagmate ', stal_ID_1{1,j}])
            figure(112);hold on
            subplot(4,5,j)
            autocorr(Data_modified_1(:,2),365)
            title(['Stalagmate ', stal_ID_1{1,j}])
        end
    end
end

% Site 2
avg_drip_2 = zeros(no_of_stalagmates2,5);
% Variance_2 = zeros(no_of_stalagmates2,1);
std_2 = zeros(no_of_stalagmates2,5);
cov_2 = zeros(no_of_stalagmates2,5);
for i = 1:5
    for j=1:no_of_stalagmates2
        Data2(:,2) = Drip2(:,1+(j-1)*2);
        %     Data2(:,3) = Drip2(:,2+(j-1)*2);
        
        %     Data_modified_2 = Data2;
        Data_modified_2 = zeros(size(Data2,1),2);
        Data_modified_2(:,1) = Data2(:,1);
        Data_modified_2(:,2) = filter( ones(k(i),1), 1, Data2(:,2) );
        Data_modified_2 = Data_modified_2(k(i):k(i):end,:);
        
        Data_modified_2(isnan(Data_modified_2(:,2)),:) = [];
        
        total_drip_2 = sum(Data_modified_2(:,2));
        total_time_2 = size(Data_modified_2,1)/4*k(i);
        % total_time = Data_modified(end,1);
        avg_drip_2(j,i) = total_drip_2/total_time_2/4;
        %     fprintf('Average drip rate per 15 minutes for %s Stalagmate = %0.1f \n', stal_ID_2{1,j}, avg_drip_2(j))
        
        %     Variance_2(j) = var(Data_modified_2(:,2));
        %     fprintf('Variance in 15 minutes drip rate for %s Stalagmate = %0.2f \n', stal_ID_2{1,j}, Variance_2(j))
        std_2(j,i) = std(Data_modified_2(:,2));
        %         cov_2(j,i) = std_2(j,i)/avg_drip_2(j,i)*100;
        cov_2(j,i) = std_2(j,i)/mean(Data_modified_2(:,2))*100;
        
        %ACF plots with daily frequency
        if i==3
            figure(14);hold on
            subplot(4,5,j)
            autocorr(Data_modified_2(:,2),200)
            title(['Stalagmate ', stal_ID_2{1,j}])
            figure(114);hold on
            subplot(4,5,j)
            autocorr(Data_modified_2(:,2),365)
            title(['Stalagmate ', stal_ID_2{1,j}])
        end
    end
end

% cov(1:size(cov_1,1),1) = cov_1;
% cov(1+size(cov_1,1):size(cov_1,1)+size(cov_2,1),1) = cov_2;
log_cov_1 = log(cov_1);
log_cov_2 = log(cov_2);

%%
figure(22);clf;hold on;
% subplot(1,2,1);hold on;
for i=1:5
    scatter(log(k(i)*15)*ones(size(avg_drip_1,1),1),log_cov_1(:,i),50,flow_class(1:size(avg_drip_1,1))','d','LineWidth',1.5);
end
for j=1:no_of_stalagmates1
    [xx,yy] = smoothLine(log(k(:)*15),log_cov_1(j,:));
    plot(xx,yy,':r','LineWidth',1.5);
%     plot(log(k(:)*15),cov_1(j,:));
end
% axis equal square
% axis tight; box on;
% xlabel('ln [Sampling frequency(min)]')
% ylabel('ln [Coefficient of Variation (COV)]')
% legend('Soda-straw','Icicle','Combined','Fracture')
% title('Sampling frequency vs Coefficient of Variation (Site 1)');

% figure(23);clf;
% subplot(1,2,2);hold on;
for i=1:5
    scatter(log(k(i)*15)*ones(size(avg_drip_2,1),1),log_cov_2(:,i),50,flow_class(size(avg_drip_1,1)+1:end)','o','LineWidth',1.5);
end
for j=1:no_of_stalagmates2
    [xx,yy] = smoothLine(log(k(:)*15),log_cov_2(j,:));
    plot(xx,yy,'--g','LineWidth',1.2);
%     plot(log(k(:)*15),cov_1(j,:));
end
axis equal square
axis tight; box on;
% xlim([0 105])
xlabel('ln [Sampling frequency(min)]')
ylabel('ln [Coefficient of Variation (COV)]')
legend('Soda-straw','Icicle','Combined','Fracture','Chamber 1','Chamber 2')
title('Sampling frequency vs Coefficient of Variation');

%% Discharge VS COV plot
figure(21);clf;hold on;
% subplot(2,2,1);
% scatter(elevation_data,avg_drip,1,'filled');
% scatter(cov_1,log(water_volume_1/365/24/3600),20,6*ones(size(water_volume_1,2),1),'filled','d');
% scatter(cov_2,log(water_volume_2/365/24/3600),20,5*ones(size(water_volume_2,2),1),'filled','o');
% scatter(cov,log(water_volume/365/24/3600),50,flow_class,'LineWidth',2);
scatter(cov_1(:,3),log(water_volume_1/365/24/3600),30,flow_class(1:size(avg_drip_1,1))','d','LineWidth',1.5);
scatter(cov_2(:,3),log(water_volume_2/365/24/3600),30,flow_class(size(avg_drip_1,1)+1:end)','o','LineWidth',1.5);
% semilogy(cov_1,water_volume_1/365/24/3600)
% semilogy(cov_2,water_volume_2/365/24/3600)
% semilogy(cov,water_volume/365/24/3600)
line([25,25], [-12.1,-13]); line([25,110], [-13,-13]);
line([0,110],[-15.37,-15.37]); line([50,50], [-13,-15.37]); line([0,50], [-13.4,-13.4]); 
line([0,25], [-12.1,-12.1]); % line([110,110], [-12.1,-10]);
% line([0,300],[-15.37,-15.37]); line([60,60], [-12.1,-15.37]); line([0,60], [-13.4,-13.4]); line([0,110], [-12.1,-12.1]); line([110,110], [-12.1,-10]);
% line([50,50], [-10,-18]); line([0,300], [-13.8,-13.8]); line([50,150], [-12.9,-12.9]); line([150,150], [-13.8,-10]);
axis equal square
axis tight; box on;
xlim([0 105]); ylim([-18 -10])
ylabel('ln [Average Discharge (l/s)]')
xlabel('Coefficient of Variation (COV)')
legend('Site 1','Site 2','Flow type')
title('ln (Discharge) vs Coefficient of Variation');
% r2_1 = (corr(log(avg_drip_1),elv1))^2;
% r2_2 = (corr(log(avg_drip_2),elv2))^2;
r5 = (corr(log(water_volume/365/24/3600),cov))^2;
r5_m = (corr(log(water_volume(flow_class==1)/365/24/3600),cov(flow_class==1)))^2;
r5_f = (corr(log(water_volume(flow_class==2)/365/24/3600),cov(flow_class==2)))^2;
r5_c = (corr(log(water_volume(flow_class==3)/365/24/3600),cov(flow_class==3)))^2;
r5_s = (corr(log(water_volume(flow_class==4)/365/24/3600),cov(flow_class==4)))^2;

% %% Friederich and Smart diagram draw
% Drip1_max = max(Drip1);
% Drip2_max = max(Drip2);
% water1_max = (Drip1_max*0.1433/1000/15/60)';
% water2_max = (Drip2_max*0.1433/1000/15/60)';

    
%%
% Site 1
active_stal_1 = 0.57;
individual_flow_count_1 = [1649 17 3 240]*2*active_stal_1; % Count of Flow type 1(M), 2(F), 3(C) and 4(S)
total_stal_count_1 = sum(individual_flow_count_1);
recharge_volume_1 = [mean(water_volume_1(1,[1:6 8 11])) mean(water_volume_1(1,10)) mean(water_volume_1(1,9)) mean(water_volume_1(1,7))]; % Assign different flow types dripping rate

individual_flow_recharge_1 = zeros(size(individual_flow_count_1));
for i = 1 : size(individual_flow_count_1,2)
    individual_flow_recharge_1(1,i) = individual_flow_count_1(1,i)*recharge_volume_1(1,i);
end
total_recharge_1 = sum(individual_flow_recharge_1);
fprintf('Site 1 Total recharge volume = %0.2f (L/yr) \n', total_recharge_1)

% site 2
active_stal_2 = 0.76;
individual_flow_count_2 = [1368 30 26 599]*2*active_stal_2; % Count of Flow type 1(M), 2(F), 3(C) and 4(S)
total_stal_count_2 = sum(individual_flow_count_2);

recharge_volume_2 = [mean(water_volume_2([1 2 6 10 12 14 15 16 18])) mean(water_volume_2([4 7 9 17])) mean(water_volume_2([3 8 11])) mean(water_volume_2([5 13]))]; % Assign different flow types dripping rate

individual_flow_recharge_2 = zeros(size(individual_flow_count_2));
for i = 1 : size(individual_flow_count_2,2)
    individual_flow_recharge_2(1,i) = individual_flow_count_2(1,i)*recharge_volume_2(1,i);
end
total_recharge_2 = sum(individual_flow_recharge_2);
fprintf('Site 2 Total recharge volume = %0.2f (L/yr) \n', total_recharge_2)

%% Take only Better quality Drip data (16 point binomial smoothing)
% Site 1
DataTable_1 = zeros(size(Drip1,1),size(Drip1,2)/2);
for j=1:size(Drip1,2)/2
    DataTable_1(:,j) = Drip1(:,2*j);
end
for j=1:size(Drip1,2)/2
    DataTable_1(isnan(DataTable_1(:,j)),:) = [];
end
% ind = find(isnan(DataTable(:,4)));

no_of_stalagmates1 = size(DataTable_1,2);
fprintf('Total Number of Stalagmates for site 1 analysis = %0.0f \n', no_of_stalagmates1)

% Site 2
DataTable_2 = zeros(size(Drip2,1),size(Drip2,2)/2);
for j=1:size(Drip2,2)/2
    DataTable_2(:,j) = Drip2(:,2*j);
end
for j=1:size(Drip2,2)/2
    DataTable_2(isnan(DataTable_2(:,j)),:) = [];
end
% ind = find(isnan(DataTable(:,4)));
no_of_stalagmates2 = size(DataTable_2,2);
fprintf('Total Number of Stalagmates for site 2 analysis = %0.0f \n', no_of_stalagmates2)


% %% Clustering
% Data12 = Data1(1:end-4,:);
% Drip12 = zeros(size(Data12,1),size(Drip1,2)/2+size(Drip2,2)/2);
% for i = 1:size(Drip1,2)/2
%     Drip12(:,i) = Drip1(1:end-4,2*i);
% end
% for i = 1:size(Drip2,2)/2
%     Drip12(:,i+size(Drip1,2)/2) = Drip2(13881:end,2*i);
% end
% 
% max_Data12 = max(Data12);
% min_Data12 = min(Data12);
% dt_Data12 = 0.0104; %equivalent of 15 minutes
% 
% %interpolate: di and vi are the interpolated dates and values
% d_Data12 = min_Data12:dt_Data12:max_Data12;
% for i=1:size(Drip12,2)
%     Drip12_i(:,i) = interp1(Data12(:,1),Drip12(:,i),d_Data12);
%     Data12_i(:,i) = d_Data12;
% end
% for j=1:size(Drip12_i,2)
%     Drip12_i(isnan(Drip12_i(:,j)),:) = [];
% end
% % DataTable = zeros(size(Drip12,1),size(Drip12,2));
% % for j=1:size(Drip12,2)
% %     DataTable(:,j) = Drip12(:,j);
% % end
% % for j=1:size(Drip12,2)
% %     DataTable(isnan(DataTable(:,j)),:) = [];
% % end
% RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
% Drip12_i = Drip12_i';
% [idx1,C1,sumd1,D1] = kmeans(Drip12_i,4);
% idx1 = idx1';
% names = stal_ID_1;
% names(end+1:size(stal_ID_1,2)+size(stal_ID_2,2)) = stal_ID_2;

%% interpolate on grid
date = cell(1,size(Drip1,2)/2+size(Drip2,2)/2);
for i=1:size(Drip1,2)/2
    date{i} = Data1;
end
for i=size(Drip1,2)/2+1:size(Drip1,2)/2+size(Drip2,2)/2
    date{i} = Data2(:,1);
end
v = cell(1,size(Drip1,2)/2+size(Drip2,2)/2);
for i=1:size(Drip1,2)/2
    v{i} = Drip1(:,2*i);
end
j=0;
for i=size(Drip1,2)/2+1:size(Drip1,2)/2+size(Drip2,2)/2
    j=j+1;
    v{i} = Drip2(:,2*j);
end
% v{1:size(Drip1,2)/2} = Drip1;
% v{size(Drip1,2)/2+1:size(Drip1,2)/2+size(Drip2,2)/2} = Drip2;

for i=1:size(v,2)
    ma(i)=max(date{i});
    mi(i)=min(date{i});
end
ma=max(ma);
mi=min(mi);

dt=0.0104; %equivalent of 15 minutes

%interpolate: di and vi are the interpolated dates and values
d=mi:dt:ma;
for i=1:size(v,2)
    vi(:,i)=interp1(date{i},v{i},d);
    di(:,i)=d;
end

% vi_r = vi(41119:end,:)';
% [idx1,C1,sumd1,D1] = kmeans(vi_r,4);

%% compute correlations and distances
for i=1:size(v,2)
    for j=1:i
        ind1=isfinite(vi(:,j));
        ind2=isfinite(vi(:,i));
        ind=ind1 & ind2;
        v1=vi(ind,j); v2=vi(ind,i);
        
        %correlation coefficients
        c=corrcoef([vi(ind,j),vi(ind,i)]);
        corr(i,j)=c(1,2);
        corr(j,i)=corr(i,j);
        
        %euclidean distance
        %v1=v1-mean(v1); v1=v1/std(v1);
        %v2=v2-mean(v2); v2=v2/std(v2);
        edist(i,j)=sqrt(mean(v1-v2).^2);
        edist(j,i)=edist(i,j);
        
        %distance in offset
        [xcf,lags,bounds] = crosscorr(v1,v2);
        m=find(xcf==max(xcf));
        offsetdist(i,j)=abs(lags(m(1)));
        offsetdist(j,i)=offsetdist(i,j);
        
        maxcorr(i,j)=max(xcf);
        maxcorr(j,i)=maxcorr(i,j);
        
%         %distance between spectrums 
%         spectrdist(i,j)=sqrt(mean(allspectr{j}(1:257,2)-allspectr{i}(1:257,2)).^2);
%         spectrdist(j,i)=spectrdist(i,j);
        
    end
    corr(i,i)=1;
    maxcorr(i,i)=1;
end
cdist=-corr+1;
maxcdist=-maxcorr+1;

%% chosing the right distance
%D=offsetdist;
% D=edist;
D=cdist;
% D=maxcdist;
% D=spectrdist;
%indtokeep=[1:22 32];
%D=D(indtokeep,indtokeep);

figure(103);clf
imagesc(D);axis square tight;colorbar

% Clustering
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));

% for i=1:size(stal_ID_1,2)-9;
%     names{i}=num2str(fileids(i));
% end
% for i=size(filename,2)-9:size(filename,2);
%     names{i}=filename{i};
% end
names = stal_ID_1;
names(end+1:size(stal_ID_1,2)+size(stal_ID_2,2)) = stal_ID_2;

% MDS
[Y,E]=cmdscale(D);

%kmeans
[IDX,C] = kmeans(Y, 4,'OnlinePhase','on','replicates',10);
IDX_R = IDX';

figure(104);clf
scatter3(Y(:,1),Y(:,1),Y(:,1),[],IDX,'o','filled')
text(Y(:,1)*1.1,Y(:,1)*1.1,Y(:,1)*1.1,names)
axis equal off

figure(105);clf
silhouette(Y,IDX)

figure(106);clf
M=linkage(squareform(D));
dendrogram(M,'labels',names)


%%
RandStream.setDefaultStream(RandStream('mt19937ar','seed',1));
cov(1:size(cov_1,1),1) = cov_1;
cov(1+size(cov_1,1):size(cov_1,1)+size(cov_2,1),1) = cov_2;

[idx1,C1,sumd1,D1] = kmeans(v,4);


% DataTable(1:size(DataTable_1,1),1) = DataTable_1;
% DataTable(1+size(DataTable_1,1):size(DataTable_1,1)+size(DataTable_2,1),1) = DataTable_2;
% 
% X1 = (DataTable_1)';
% X2 = (DataTable_2)';
% [idx1,C1,sumd1,D1] = kmeans(X1,4);
% idx1 = idx1'; idx2 = idx2';
% [idx2,C2,sumd2,D2] = kmeans(X2,4);
% % figure(71);clf;
% % plot(X(:,1),X(:,2),'k*','MarkerSize',5);
% % title 'Fisher''s Iris Data';
% % xlabel 'Petal Lengths (cm)';
% % ylabel 'Petal Widths (cm)';

% T = clusterdata(DataTable_1(1:100,1:3),2);
% figure(51);clf;
% scatter3(DataTable_1(1:100,1),DataTable_1(1:100,2),DataTable_1(1:100,3),30,T,'filled')

%% Plot Correlation Matrix with Elevation
% elevation = [22.39 22.36 22.35 22.46 22.54 21.85 22.53 21.90 22.08 21.80 21.83];
% Site 1
[elv_ind_1 I1] = sort(flow_class(1:size(avg_drip_1,1)));
stal_ID_elv_1 = stal_ID_1(I1);
DataTable_elv_1 = DataTable_1(:,I1);

[r1,~] = corr(DataTable_elv_1,'type','Kendall');
% r = kendall(DataTable);
figure(22);clf;
% pcolor(r);
imagesc(r1);
n = size(r1,1);
axis equal square
axis tight
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', stal_ID_elv_1); % set x-axis labels
set(gca, 'YTickLabel', stal_ID_elv_1); % set y-axis labels
title('Kendall Tau Correlation matrix of Site 1 drip data'); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1]); % set the scale of colorbar

% Site 2
[elv_ind_2 I2] = sort(flow_class(size(avg_drip_1,1)+1:end));
stal_ID_elv_2 = stal_ID_2(I2);
DataTable_elv_2 = DataTable_2(:,I2);

[r2,~] = corr(DataTable_elv_2,'type','Kendall');
% r = kendall(DataTable);
figure(23);clf;
% pcolor(r);
imagesc(r2);
n = size(r2,1);
axis equal square
axis tight
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', stal_ID_elv_2); % set x-axis labels
set(gca, 'YTickLabel', stal_ID_elv_2); % set y-axis labels
title('Kendall Tau Correlation matrix of Site 2 drip data'); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1]); % set the scale of colorbar

%% Plot Correlation Matrix with Elevation
% elevation = [22.39 22.36 22.35 22.46 22.54 21.85 22.53 21.90 22.08 21.80 21.83];
% Site 1
[elv_ind_1 I1] = sort(depth1);
stal_ID_elv_1 = stal_ID_1(I1);
DataTable_elv_1 = DataTable_1(:,I1);

[r1,~] = corr(DataTable_elv_1,'type','Kendall');
% r = kendall(DataTable);
figure(22);clf;
% pcolor(r);
imagesc(r1);
n = size(r1,1);
axis equal square
axis tight
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', stal_ID_elv_1); % set x-axis labels
set(gca, 'YTickLabel', stal_ID_elv_1); % set y-axis labels
title('Kendall Tau Correlation matrix of Site 1 drip data'); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1]); % set the scale of colorbar

% Site 2
[elv_ind_2 I2] = sort(depth2);
stal_ID_elv_2 = stal_ID_2(I2);
DataTable_elv_2 = DataTable_2(:,I2);

[r2,~] = corr(DataTable_elv_2,'type','Kendall');
% r = kendall(DataTable);
figure(23);clf;
% pcolor(r);
imagesc(r2);
n = size(r2,1);
axis equal square
axis tight
set(gca, 'XTick', 1:n); % center x-axis ticks on bins
set(gca, 'YTick', 1:n); % center y-axis ticks on bins
set(gca, 'XTickLabel', stal_ID_elv_2); % set x-axis labels
set(gca, 'YTickLabel', stal_ID_elv_2); % set y-axis labels
title('Kendall Tau Correlation matrix of Site 2 drip data'); % set title
colormap('jet'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1]); % set the scale of colorbar

%% Plot Correlation coefficients
% figure(105);clf;
% corrplot(DataTable_2(:,[10 11 12]),{'1vii' '1viii' '1x'},'a')
% ylim([-1,1])
% title('Correlation matrix between fracture flow drips time series');
% 
% figure(106);clf;
% corrplot(DataTable_2(:,[1:6 8]),{'1A' '1B' '1i' '1ii' '1iii' '1iv' '1v'},'a')
% title('Correlation matrix between Site 1 left patch drips time series');


