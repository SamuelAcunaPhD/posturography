% Filename: demo_compare.m
% Author:   Samuel Acuña
% Date:     28 Feb 2019 
% Description:
% Developed for kids, quickly visualize and test two collections
% 
% collect only: TIME, COPx, COPy


clear; close all; clc;

[filename, pathname] = uigetfile('*.csv', 'pick two files','MultiSelect', 'on');
if all(class(filename)=='char')
    nTrials = 1;
    filename = {filename};
else
    nTrials = length(filename);
end

% convert units
conversion = 100; %cm
units = 'cm';

XMIN = 0; XMAX = 0; YMIN = 0; YMAX = 0;
for tr = 1:nTrials % cycle through every trial

    
    file = [pathname filename{tr}];
    
    
    % extract data from balance plate. columns are: seq time Fz Mx My CoPx CoPy
    M = csvread(file,1);
    
    % assign columns to variables
    time  = M(:,1);
    COPml1 = conversion*detrend(M(:,2)); 
    COPap1 = conversion*-detrend(M(:,3)); % need to invert, since forward is negative.
        
    % filter everything above 10 Hz
    sampRate = 1/(time(2)-time(1)); % find sample rate
    cutoffFrequency = 10;
    [b,a]=butter(4,cutoffFrequency/(sampRate/2)); %butterworth filter, 4th order
    COPml2{tr} = filtfilt(b,a,COPml1);
    COPap2{tr} = filtfilt(b,a,COPap1);
    
    totalTime{tr} = max(round(time,1));
    
    % store limits for plotting
    if min(COPml1) < XMIN
        XMIN = min(COPml1);
    end
    if min(COPap1) < YMIN
        YMIN = min(COPap1);
    end
    if max(COPml1) > XMAX
        XMAX = max(COPml1);
    end
    if max(COPap1) > YMAX
        YMAX = max(COPap1);
    end
    
end

figure();
for tr = 1:nTrials % cycle through every trial
    COPml = COPml2{tr};
    COPap = COPap2{tr};
    
    % 1: plot COP
    ax1 = subplot(nTrials,3,3*tr-2);
    plot(COPml,COPap,'.');
    hold on
    plot(COPml(1),COPap(1),'gs');
    plot(COPml(end),COPap(end),'rs');
    hold off
    xlabel(['Side to Side [' units ']']);
    ylabel(['Front to back [' units ']']);
    axis equal
    axis(1.6*[XMIN XMAX YMIN YMAX]);
    title('Center of Pressure');
    axis square
    ax1.FontSize = 12;
    text(mean([XMIN XMAX]),1.35*YMAX,filename{tr}(1:end-4),'Color','red','FontSize',12,'HorizontalAlignment','center','Interpreter','none');
    text(mean([XMIN XMAX]),1.35*YMIN,[ num2str(totalTime{tr}) ' seconds'],'Color','red','FontSize',12,'HorizontalAlignment','center');
    
    % 2: plot area
    d = [COPml, COPap];       %zero mean COP values
    [m,n]=size(d);         %returns m=rows, n=columns of d
    mean_d=mean(d);        %returns row vector with column means of d
    cov_mtx=cov(d);        %covariance matrix for d
    [V,D]=eig(cov_mtx);    %V=eigenvectors, D=eigenvalues of cov_mtx
    semimaj=[mean_d; mean_d+2.45*sqrt(D(1,1))*V(:,1)']; %center and end of semimajor axis
    semimin=[mean_d; mean_d+2.45*sqrt(D(2,2))*V(:,2)']; %center and end of semiminor axis
    theta=linspace(0,2*pi,41)';
    ellipse=2.45*sqrt(D(1,1))*cos(theta)*V(:,1)' + 2.45*sqrt(D(2,2))*sin(theta)*V(:,2)' + ones(size(theta))*mean_d;
    c = 5.991; % for 95%. For 90%, use c = 4.605
    Area95 = c*pi*sqrt(D(1,1)*D(2,2));
    
    ax2 = subplot(nTrials,3,3*tr-1)
    plot(d(:,1),d(:,2),'.');  %scatter plot with x=column 1 of d, y=column 2
    hold on;
    plot(semimaj(:,1),semimaj(:,2),'r','LineWidth',2);
    plot(semimin(:,1),semimin(:,2) ,'r','LineWidth',2);
    plot(ellipse(:,1),ellipse(:,2) ,'g','LineWidth',2);
    axis equal
    axis(1.6*[XMIN XMAX YMIN YMAX]);
    axis square
    text(mean([XMAX XMIN]),1.35*YMAX,[num2str(round(Area95,1)) ' ' units],'Color','red','FontSize',12,'HorizontalAlignment','center');
    title('Area');
    ax2.FontSize = 12;
    
    % 3: plot range
    ax3 = subplot(nTrials,3,3*tr)
    plot(COPml,COPap,'.');
    rangeCOPml = max(COPml)-min(COPml);
    rangeCOPap = max(COPap)-min(COPap);
    % ML line
    line([min(COPml) max(COPml)], 1.1*[YMIN YMIN],'LineWidth',2,'Color','r');
    line([min(COPml) min(COPml)], [1.1*YMIN YMIN],'LineWidth',2,'Color','r');
    line([max(COPml) max(COPml)], [1.1*YMIN YMIN],'LineWidth',2,'Color','r');
    % AP line
    line(1.1*[XMIN XMIN], [min(COPap) max(COPap)],'LineWidth',2,'Color','r');
    line([1.1*XMIN XMIN],[min(COPap) min(COPap)],'LineWidth',2,'Color','r');
    line([1.1*XMIN XMIN],[max(COPap) max(COPap)],'LineWidth',2,'Color','r');
    % text
    text(0,1.35*YMIN,[num2str(round(rangeCOPml,1)) ' ' units],'Color','red','FontSize',12,'HorizontalAlignment','center');
    text(1.35*XMIN,0,[num2str(round(rangeCOPap,1)) ' ' units],'Color','red','FontSize',12,'HorizontalAlignment','center','Rotation',90);
    axis equal
    axis(1.6*[XMIN XMAX YMIN YMAX]);
    axis square
    title('Range');
    ax3.FontSize = 12;
    
end