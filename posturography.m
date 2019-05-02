% Filename: posturography.m
% Author:   Samuel Acuña
% Updated:  27 Nov 2018
% Description:
% functions used to compare posturography. This code should serve as a
% starting point for all your posturography needs. Note this is a static
% class, so you dont need a constructor to access it.
%
% note: 06 Mar 2019
% I have changed 'detrend' in COP functions to just remove the DC mean, not a linear trend. Not sure if this should be changed back
%
% run this code using a driver file. For example:
%     filenames = {'eyesOpen_1.csv';'eyesOpen_2.csv';'eyesOpen_3.csv';'eyesClosed_1.csv';'eyesClosed_2.csv';'eyesClosed_3.csv'};
%     COP = posturography.load(filenames); % loads data files, filters, detrends
%     posturography.statokinesiogram_compare(COP,[1 2 3; 4 5 6]); % plot multiple statokinesiograms
%     posturography.stabilogram(COP,[1 2 3; 4 5 6]); % plot multiple stabilograms in specified order
%     posturography.estimateCOG(COP,[1 2 3; 4 5 6]); % plot Center of Gravity
%     [~,metrics] = posturography.calcMetrics(COP(1)) % just output metrics
%     COP = posturography.calcMetrics(COP); % calculate common metrics
%     posturography.plotMetrics(COP,[1 2 3; 4 5 6]);
%     posturography.plotSpectralAnalysis(COP(1)) % single FFT plot
%     posturography.plotSpectralAnalysis(COP,[1 2 3; 4 5 6]) % all FFTs

classdef posturography
    properties (Constant, Access = 'public')
        % nothing yet
    end
    methods (Static)
        function COP = load(filenames)
            % loads posturography raw data
            for i = 1:length(filenames)
                
                %note: some force plates will require you to calculate COP
                %manually.
                
                % extract data from balance plate.
                M = csvread(filenames{i},1);
                
                % make sure I have correct amount of columns (seven, in this case)
                % columns are: seq time Fz Mx My CoPx CoPy
                [nRows, nCols] = size(M);
                assert(nCols == 7);
                
                % assign columns to variables
                %seq  = M(:,1);
                time = M(:,2);
                Fz   = M(:,3); % force in Z direction
                %Mx   = M(:,4); % moment about x axis
                %My   = M(:,5); % moment about y axis
                %COPx = detrend(M(:,6)); % center of pressure x direction (ML direction usually)
                %COPy = detrend(M(:,7)); % center of pressure y direction (AP direction usually)
                COPx = detrend(M(:,6),'constant'); % center of pressure x direction (ML direction usually)
                COPy = detrend(M(:,7),'constant'); % center of pressure y direction (AP direction usually)
                
                sampRate = 1/(time(2)-time(1)); % find sample rate
                mass = mean(Fz)/9.81; % newtons to kg
                
                %% filter data
                % filter everything above 10 Hz, recommended by 1. Ruhe A, Fejer R, Walker B. The test-retest reliability of centre of pressure measures in bipedal static task conditions - A systematic review of the literature. Gait Posture 32: 436–445, 2010.
                cutoffFrequency = 10;
                [b,a]=butter(4,cutoffFrequency/(sampRate/2)); %butterworth filter, 4th order
                COPx = filtfilt(b,a,COPx);
                COPy = filtfilt(b,a,COPy);
                
                %% zero mean COP
                % optional, might not want to do when comparing different conditions on
                % same statokinesiogram
                
                %COPx = detrend(COPx);
                %COPy = detrend(COPy);
                COPx = detrend(COPx,'constant');
                COPy = detrend(COPy,'constant');
                
                %% retain extracted data
                COP(i).filename = filenames{i};
                COP(i).time = time;
                COP(i).Fz = Fz;
                COP(i).ml = COPx;
                COP(i).ap = COPy;
                COP(i).mass = mass;
                COP(i).sampRate = sampRate;
                
                
            end
        end
        function COP = load_twoforceplates(data)
            % for loading posturography data that is
            % returns COP data for each foot and net COP
            % data is filtered and shifted to have mean(netCOP) = 0,0
            % COP output as a structure
            
            % extract raw data
            FZ1 = data.FP_FZ1;
            FZ2 = data.FP_FZ2;
            COPx1 = data.FP_COPx1; % left foot
            COPx2 = data.FP_COPx2; % right foot
            COPy1 = data.FP_COPy1; 
            COPy2 = data.FP_COPy2;
            
            % find net COP
            % D.A. Winter, F. Prince, P. Stergiou, C. Powell, Medial-Lateral and Anterior-Posterior Motor-Responses Associated with Center of Pressure Changes in Quiet Standing, Neurosci. Res. Commun. 12 (1993) 141?148.
            FZnet = FZ1+FZ2;
            FZ1_ratio = FZ1./FZnet;
            FZ2_ratio = FZ2./FZnet;
            COPml = COPx1.*FZ1_ratio+COPx2.*FZ2_ratio;
            COPap = COPy1.*FZ1_ratio+COPy2.*FZ2_ratio;
            
            if any(isnan(COPap))
                disp('COPap: NaNs detected. interpolating over them');
                COPap = inpaint_nans(COPap);
            end
            if any(isnan(COPml))
                disp('COPml: NaNs detected. interpolating over them');
                COPml = inpaint_nans(COPml);
            end
            
            sampRate = data.fs; % sample rate
            
            % filter everything above 10 Hz, recommended by 1. Ruhe A, Fejer R, Walker B. The test-retest reliability of centre of pressure measures in bipedal static task conditions - A systematic review of the literature. Gait Posture 32: 436–445, 2010.
            cutoffFrequency = 10;
            [b,a]=butter(4,cutoffFrequency/(sampRate/2)); %butterworth filter, 4th order
            COPml = filtfilt(b,a,COPml);
            COPap = filtfilt(b,a,COPap);
            
            % find distances to have mean COP at zero
            %SHIFTml = COPml - detrend(COPml);
            %SHIFTap = COPap - detrend(COPap);
            SHIFTml = COPml - detrend(COPml,'constant');
            SHIFTap = COPap - detrend(COPap,'constant');
            
            % detrend COP data
            %COPml = detrend(COPml);
            %COPap = detrend(COPap);
            COPml = detrend(COPml,'constant');
            COPap = detrend(COPap,'constant');
            
            % filter components of COP
            COPx1 = filtfilt(b,a,COPx1); COPx2 = filtfilt(b,a,COPx2);
            COPy1 = filtfilt(b,a,COPy1); COPy2 = filtfilt(b,a,COPy2);
            
            % shift data to have mean COP at zero
            COPx1 = COPx1 - SHIFTml; COPx2 = COPx2 - SHIFTml;
            COPy1 = COPy1 - SHIFTap; COPy2 = COPy2 - SHIFTap;
            
            % output
            COP.ml = COPml;
            COP.ap = COPap;
            COP.x1 = COPx1; % left foot
            COP.x2 = COPx2; % right foot
            COP.y1 = COPy1;
            COP.y2 = COPy2;
            COP.FZ1_ratio = FZ1_ratio;
            COP.FZ2_ratio = FZ2_ratio;
        end
        function fig = statokinesiogram(COP, order)
            % INPUT: 
            %      COP: structure, such that COP.ml and COP.ap exist
            %      order: Optional. array specifying how to plot elements of COP
            % 
            % for example: if length(COP) = 6; order = [1 2 3; 4 5 6]; then
            % subplot will be 2x3
            
            % default values
            if nargin == 1
                order = 1:length(COP);
            end
            
            assert(numel(order)==length(COP)); 
            
            % setup figures
            fig(1) = figure();
            
            % plotting options
            plot_size = size(order);
            order2 = reshape(order,1,numel(order));
            lineWidth = 2;
            fontSize = 12;
            unitConversion=1000;
            units = 'mm';
            
            % set axis size 
            axis_size = [0 0 0 0];
            for i = 1:length(COP)
                if min(COP(i).ml) < axis_size(1)
                    axis_size = [min(COP(i).ml) axis_size(2:4)];
                end
                if max(COP(i).ml) > axis_size(2)
                    axis_size = [axis_size(1) max(COP(i).ml) axis_size(3:4)];
                end
                if min(COP(i).ap) < axis_size(3)
                    axis_size = [axis_size(1:2) min(COP(i).ap) axis_size(4)];
                end
                if max(COP(i).ap) > axis_size(4)
                    axis_size = [axis_size(1:3) max(COP(i).ap)];
                end
            end
            
            
            
            
            for i = order2
                ax(i) = subplot(plot_size(1),plot_size(2),i)
                
                % convert units
                COPml = COP(i).ml*unitConversion;
                COPap = COP(i).ap*unitConversion;
                
                % plot
                plot(COPml,COPap,'LineWidth',lineWidth)
                if all(isfield(COP,{'x1','y1','x2','y2'}))
                    COPx1 = COP(i).x1*unitConversion;
                    COPy1 = COP(i).y1*unitConversion;
                    COPx2 = COP(i).x2*unitConversion;
                    COPy2 = COP(i).y2*unitConversion;
                    hold on
                    plot(COPx1,COPy1,'LineWidth',lineWidth)
                    plot(COPx2,COPy2,'LineWidth',lineWidth)
                    hold off
                    legend('Net COP','Left','Right')
                end
                
                % axes
                axis(axis_size*unitConversion);
                axis equal
                
                % labels
                xlabel(['COP_M_L (' units ')'],'FontSize',fontSize)
                ylabel(['COP_A_P (' units ')'],'FontSize',fontSize)
                
                % title
                if isfield(COP,'filename')
                    title(['COP: ' COP(i).filename],'Interpreter','none','FontSize',fontSize);
                else
                    title('Center of Pressure','FontSize',fontSize);
                end
            end
            linkaxes(ax,'xy'); % link axes so they zoom together
            
            % save figure info
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8.5 11];
            fig.Name = ['statokinesiogram'];

        end
        function statokinesiogram_animation(COP)
            % animate single COP trajectory as statokinesiogram
            % INPUT:
            %      COP: structure, such that COP.ml and COP.ap exist
            
            if length(COP)>1
                error('statokinesiogram_animation requires only one COP structure.')
            end
            % setup figures
            fig(1) = figure();
            
            % set plotting options
            axis_size = [min(min([COP.ml])),max(max([COP.ml])),min(min([COP.ap])),max(max([COP.ap]))];
            lineWidth = 2;
            fontSize = 12;
            unitConversion=1000;
            units = 'mm';
            
            % axes
            ax = axes;
            axis(axis_size*unitConversion);
            axis equal
            xlimits = ax.XLim;
            ax.XLim = 1.1*xlimits; % set limits to 110%
            
            % labels
            xlabel(['COP_M_L (' units ')'],'FontSize',fontSize)
            ylabel(['COP_A_P (' units ')'],'FontSize',fontSize)
            
            % title
            if isfield(COP,'filename')
                title(['COP: ' COP.filename],'Interpreter','none','FontSize',fontSize);
            else
                title('Center of Pressure','FontSize',fontSize);
            end
            
            
            % convert units
            COPml = COP.ml*unitConversion;
            COPap = COP.ap*unitConversion;
            
            
            % plot animation
            hold(ax)
            cometLength = 0.1;
            comet(COPml,COPap,cometLength)
            
            % save figure info
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8.5 11];
            fig.Name = ['statokinesiogram'];
        end
        function fig = stabilogram(COP,order)
            % plots stabilgram (COP over time)
            % INPUT: 
            %      COP: structure, such that COP.ml,COP.ap, COP.time exist
            %      order: Optional. array specifying how to plot elements of COP
            % 
            % for example: if length(COP) = 6; order = [1 2 3; 4 5 6]; then
            % subplot will be 2x3
            
            % default values
            if nargin == 1
                order = 1:length(COP);
            end
            
            assert(numel(order)==length(COP)); 
            
            % setup figures
            fig(1) = figure();
            
            % set plotting options
            ylimSize = [min([min([COP.ml]) min([COP.ap])]),max([max([COP.ml]) max([COP.ap])])];
            plot_size = size(order);
            order2 = reshape(order,1,numel(order));
            lineWidth = 2;
            fontSize = 12;
            unitConversion=1000;
            units = 'mm';
            
            for i = order2
                ax(i) = subplot(plot_size(1),plot_size(2),i)
                
                % convert units
                COPml = COP(i).ml*unitConversion;
                COPap = COP(i).ap*unitConversion;
                
                % plot
                plot(COP(i).time,COPml,'LineWidth',lineWidth)
                hold on
                plot(COP(i).time,COPap,'LineWidth',lineWidth)
                hold off
               
                
                % axes
                ylim(ylimSize*unitConversion);
                
                % labels
                xlabel(['Time (s)'],'FontSize',fontSize)
                ylabel(['COP Displacement (' units ')'],'FontSize',fontSize)
                
                % title
                if isfield(COP,'filename')
                    title(['COP: ' COP(i).filename],'Interpreter','none','FontSize',fontSize);
                else
                    title('Center of Pressure','FontSize',fontSize);
                end
            end
            
            % legend
            subplot(plot_size(1),plot_size(2),i);
            legend({'COP_M_L','COP_A_P'}, 'Location','best','boxoff')
            
            linkaxes(ax,'xy'); % link axes so they zoom together
            
            % save figure info
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8.5 11];
            fig.Name = ['stabilogram'];
        end
        function fig = stabilogram_animation(COP)
            % animate single COP trajectory as stabilogram
            % INPUT:
            %      COP: structure, such that COP.ml and COP.ap exist
            
            
            if length(COP)>1
                error('stabiliogram_animation requires only one COP structure.')
            end
            % setup figures
            fig(1) = figure();
            
            % set plotting options
            ylimSize = [min([min([COP.ml]) min([COP.ap])]),max([max([COP.ml]) max([COP.ap])])];
            xlimSize = [COP.time(1) COP.time(end)];
            lineWidth = 2;
            fontSize = 12;
            unitConversion=1000;
            units = 'mm';
            
            % axes
            ylim(ylimSize*unitConversion);
            xlim(xlimSize);
            
            
            % labels
            xlabel(['Time (s)'],'FontSize',fontSize)
            ylabel(['COP Displacement (' units ')'],'FontSize',fontSize)
            
            % title
            if isfield(COP,'filename')
                title(['COP: ' COP.filename],'Interpreter','none','FontSize',fontSize);
            else
                title('Center of Pressure','FontSize',fontSize);
            end
            
            
            % convert units
            COPml = COP.ml*unitConversion;
            COPap = COP.ap*unitConversion;
            
            
            % plot animation
            timeStep = 50;
            h(1) = animatedline; h(1).Color = [0.2539    0.4102    0.8789]; h(1).LineWidth = lineWidth;
            h(2) = animatedline; h(2).Color = [0.8594    0.0781    0.2344]; h(2).LineWidth = lineWidth;
            addpoints(h(1),COP.time(1),COPml(1));
            addpoints(h(2),COP.time(1),COPap(1));
            drawnow
            legend(h,{'COP_M_L','COP_A_P'}, 'Location','northeast','boxoff')
            for k = 1:timeStep:length(COPml)-timeStep
                addpoints(h(1),COP.time(k:(k+timeStep-1)),COPml(k:(k+timeStep-1)));
                addpoints(h(2),COP.time(k:(k+timeStep-1)),COPap(k:(k+timeStep-1)));
                drawnow
            end
            
            % save figure info
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8.5 11];
            fig.Name = ['stabilogram'];
        end
        function [fig,COG] = estimateCOG(COP,order)
            % estimate and plot center of gravity in AP direction
            % INPUT: 
            %      COP: structure, such that COP.ml,COP.ap, COP.time exist
            %      order: Optional. array specifying how to plot elements of COP
            % 
            % for example: if length(COP) = 6; order = [1 2 3; 4 5 6]; then
            % subplot will be 2x3
            %
            % estimate Center of Gravity (COG), using filtering method and
            % considering the body as an inverted pendulum. Simple and
            % fast, but not super accurate.
            % using method from:?Duarte, M., & Freitas, S. M. S. F. (2010). Revision of posturography based on force plate for balance evaluation. Revista Brasileira de Fisioterapia (São Carlos (São Paulo, Brazil)), 14(3), 183?192. http://doi.org/10.1590/S1413-35552010000300003
            
            % default values
            if nargin == 1
                order = 1:length(COP);
            end
            
            assert(numel(order)==length(COP)); 
            
            % setup figures
            fig(1) = figure();
            
            % set plotting options
            ylimSize = [min([min([COP.ap])]),max([max([COP.ap])])];
            plot_size = size(order);
            order2 = reshape(order,1,numel(order));
            lineWidth = 2;
            fontSize = 12;
            unitConversion=1000;
            units = 'mm';
            
            for i = order2
                ax(i) = subplot(plot_size(1),plot_size(2),i)
                
                
                % estimate Center of Gravity (COG), using filtering method
                cutoffFrequency = 0.5;
                [b,a]=butter(4,cutoffFrequency/(COP(i).sampRate/2)); %butterworth filter, 4th order
                COG(i).ap = filtfilt(b,a,COP(i).ap);
                
                
                % convert units
                COPap = COP(i).ap*unitConversion;
                COGap = COG(i).ap*unitConversion;
                
                % plot
                plot(COP(i).time,COPap,'LineWidth',lineWidth)
                hold on
                plot(COP(i).time,COGap,'LineWidth',lineWidth)
                hold off
               
                
                % axes
                ylim(ylimSize*unitConversion);
                
                % labels
                xlabel(['Time (s)'],'FontSize',fontSize)
                ylabel(['COP Displacement (' units ')'],'FontSize',fontSize)
                
                % title
                if isfield(COP,'filename')
                    title(['COG: ' COP(i).filename],'Interpreter','none','FontSize',fontSize);
                else
                    title('Center of Gravity','FontSize',fontSize);
                end
            end
            
            % legend
            subplot(plot_size(1),plot_size(2),i);
            legend({'COP_A_P','COG_A_P'}, 'Location','best','boxoff')
            
            linkaxes(ax,'xy'); % link axes so they zoom together
            
            % save figure info
            fig.PaperUnits = 'inches';
            fig.PaperPosition = [0 0 8.5 11];
            fig.Name = ['centerOfGravity'];
        end
        function [COP, metrics] = calcMetrics(COP)
            % wrapper for posturographyMetrics. Look at this function for
            % more detail, which I highly recommend.
            for i = 1:length(COP)
                metrics(i) = posturography.posturographyMetrics(COP(i).ml,COP(i).ap,COP(i).sampRate);
                COP(i).metrics = metrics(i);
            end
        end
        function metrics = posturographyMetrics(COPml, COPap, sampRate)
            % INPUTS:
            % COPml = center of pressure over time, mediolateral direction
            % COPap = center of pressure over time, Antero-posterior direction
            % sampRate = sampling rate of COP traces
            %
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Posturography Measurements:
            % RDavg             Average Radial Displacement, aka Mean COP distance
            % Vavg              Average Velocity
            % Vml               Average Velocity in the mediolateral direction
            % Vap               Average Velocity in the Antero-posterior direction
            % SDml              Standard deviation of the COP in the mediolateral direction
            % SDap              Standard deviation of the COP in the Antero-posterior direction
            % SDrd              Standard deviation of the Radial displacement
            % Area95            Area of the 95% confidence ellipse (with plot)
            % Area95_perSec     Area of the 95% confidence ellipse per second
            % Area95_perSec_std Area of the 95% confidence ellipse per second, standard deviation
            % Hull              Hull of COP (NOT IMPLEMENTED YET)
            % Circ95            95% prediction circumference area
            % Circ95_perSec     95% prediction circumference area per second
            % Circ95_perSec_std 95% prediction circumference area per second, standard deviation
            % TURNi             Turns index
            % Beta              mean angle deviance from AP
            % Beta_std          mean angle deviance from AP
            % Rml               Range of COP ml displacement
            % Rap               Range of COP ap displacement
            % SPml              Sway Path ml
            % SPap              Sway Path ap
            % SPr               Resultant Sway Path, very similar to RDavg. aka total radial displacement
            % PL                Pathlength, proportional to Vavg
            % PLml              Pathlength, ML
            % PLap              Pathlength, AP
            % PLn               Normalized Pathlength, proportional to Vavg
            % PLmln             Normalized Pathlength, ML
            % PLapn             Normalized Pathlength, AP
            % Fpeakml           Peak, Power spectral density ML
            % Fmeanml           Mean, Power spectral density ML
            % F50ml             Median, Power spectral density ML
            % F80ml             Frequency band that contains up to 80% of the spectrum, Power spectral density ML
            % Fpeakap           Peak, Power spectral density AP
            % Fmeanap           Mean, Power spectral density AP
            % F50ap             Median, Power spectral density AP
            % F80ap             Frequency band that contains up to 80% of the spectrum, Power spectral density AP
            % FFTfrequencies    Frequency vector from FFT
            % FFTml             ML FFT coefficients
            % FFTap             AP FFT coefficients
            
            % NOT USED:
            % RMSml, RMSap      Root Mean Square. With zero mean data, RMS is the same as SDml and SDap
            
            
            % citations:
            % 1. appendix of Baig, Dansereau, Chan, et al. 2012. Internation
            % journal of electrical and computer engineering. "Cluster analysis of
            % Center of Pressure Measures"
            % 2. Duarte, Freitas, Zatsiorsky. 2011. 'Control of Equilibrium in
            % Humans - Sway over sway'. chapter 10 of book "Motor Control:
            % Theories, Experiments, and Applications" edited by Latash and Danion
            % 3. Schubert, P., Kirchner, M., Schmidtbleicher, D., & Haas, C. T.
            % (2012). About the structure of posturography: Sampling duration,
            % parametrization, focus of attention (part I). Journal of Biomedical
            % Science and Engineering, 05(09), 496–507.
            % http://doi.org/10.4236/jbise.2012.59062
            % 4. Hufschmidt, A., Dichgans, J., Mauritz, K.-H., & Hufschmidt, M.
            % (1980). Some methods and parameters of body sway quantification and
            % their neurological applications. Archiv Für Psychiatrie Und
            % Nervenkrankheiten, 228(2), 135–150. http://doi.org/10.1007/BF00365601
            % 5. Whitney 2011. Gait & Posture. "A comparison of accelerometry and
            % center of pressure measures during computerized dynamic
            % posturography: A measure of balance"
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % check COP are same length
            if length(COPap) ~= length(COPml)
                error('COPml and COPap are different lengths');
            end
            
            % OUTPUT STRUCTURE
            metrics = struct();
            
            %% RDavg             Average Radial Displacement, aka Mean COP distance
            % note: this is very similar to the resultant sway path (SPr). dont  use both of them.
            rad = sqrt(COPml.^2 + COPap.^2);
            metrics.RDavg = mean(rad);
            
            %% Vavg              Average Velocity
            metrics.Vavg =  sum(sqrt(diff(COPap).^2+diff(COPml).^2))*sampRate/length(COPap);
            
            %% Vml               Average Velocity in the mediolateral direction
            metrics.Vml = sum(abs(diff(COPml)))*sampRate/length(COPml);
            
            %% Vap               Average Velocity in the Antero-posterior direction
            metrics.Vap = sum(abs(diff(COPap)))*sampRate/length(COPap);
            
            %% SDml              Standard deviation of the COP in the mediolateral direction
            metrics.SDml = std(COPml);
            
            %% SDap              Standard deviation of the COP in the Antero-posterior direction
            metrics.SDap = std(COPap);
            
            %% SDrd              Standard deviation of the Radial displacement
            metrics.SDrd = std(rad);
            
            %% Area95            Area of the 95% confidence ellipse
            % from the table of F statistics at a confidence level of 1-alpha
            % with alpha=0.05 when the sample size is >120.
            % F = 3.00;
            
            % sigma_ml = SDml
            % sigma_ap = SDap
            
            % correlation coefficient between COPml and COPap
            % sigma_ml_ap = 1/N*sum(COPml.*COPap)/(sigma_ml*sigma_ap)
            % sigma_corr_coeffs = corrcoef(COPml, COPap);
            % sigma_ml_ap = sigma_corr_coeffs(1,2);
            
            % area computed as Schmit et al, 2006 ? This isn't true, I dont
            % think
            % AREA_Baig(subj,i) = (2*pi*F)*sqrt(sigma_ml^2*sigma_ap^2 - sigma_ml_ap^2)
            
            % 95 confidence elipse (adapted from Frank Borg, Kokkola, Finland via University of Delaware)
            % http://www.udel.edu/biology/rosewc/kaap686/reserve/cop/center%20of%20pressure.html
            d = [COPml, COPap];       %zero mean COP values
            cov_mtx=cov(d);           %covariance matrix for COP
            [V,D]=eig(cov_mtx);       %V=eigenvectors, D=eigenvalues of cov_mtx
            
            % confidence from table of critical values for chi squared
            % distribution with 2 Degrees of Freedom
            c = 5.991; % for 95%. For 90%, use c = 4.605
            metrics.Area95 = c*pi*sqrt(D(1,1)*D(2,2));
            
            % alternate:
            % [vec, val] = eig(cov(COPap, COPml));
            % Area952 = pi*prod(2.4478*sqrt(svd(val)));
            
            % ALSO...
            % mean_d = mean(d);
            % semimaj=[mean_d; mean_d+sqrt(c*D(1,1))*V(:,1)'];
            %                           %center and end of semimajor axis
            % semimin=[mean_d; mean_d+sqrt(c*D(2,2))*V(:,2)'];
            %                          %center and end of semiminor axis
            % L_maj = sqrt(c*D(1,1)); % OR pdist(semimaj,'euclidean');
            % L_min = sqrt(c*D(2,2)); % OR pdist(semimin,'euclidean');
            % Area = pi*L_maj*L_min
            
            % points on the ellipse, just for plotting sake
            % theta=linspace(0,2*pi,41)';
            % ellipse=sqrt(c*D(1,1))*cos(theta)*V(:,1)' + sqrt(c*D(2,2))*sin(theta)*V(:,2)' + ones(size(theta))*d
            
            % % plot ellipse:
            % if (0) % set true to plot
            %     figure;
            %     if i == 1
            %         plot(d(:,1),d(:,2),'b.');  %scatter plot with x=column 1 of d, y=column 2
            %     else
            %         plot(d(:,1),d(:,2),'k.');  %scatter plot with x=column 1 of d, y=column 2
            %     end
            %
            %     hold on;
            %     plot(semimaj(:,1),semimaj(:,2),'r','LineWidth',2);
            %     plot(semimin(:,1),semimin(:,2) ,'r','LineWidth',2);
            %     plot(ellipse(:,1),ellipse(:,2) ,'g','LineWidth',2);
            %     title(filenames{i})
            %     axis square
            %     % axis([-0.08 0.015 -0.15 -0.08])
            % end
            
            % I'm getting issues with the above code for plotting. try this:
            % d = [COPml, COPap];       %zero mean COP values
            % [m,n]=size(d)         %returns m=rows, n=columns of d
            % mean_d=mean(d)        %returns row vector with column means of d
            % cov_mtx=cov(d)        %covariance matrix for d
            % [V,D]=eig(cov_mtx)    %V=eigenvectors, D=eigenvalues of cov_mtx
            % semimaj=[mean_d; mean_d+2.45*sqrt(D(1,1))*V(:,1)']
            % %center and end of semimajor axis
            % semimin=[mean_d; mean_d+2.45*sqrt(D(2,2))*V(:,2)']
            % %center and end of semiminor axis
            % theta=linspace(0,2*pi,41)';
            % ellipse=2.45*sqrt(D(1,1))*cos(theta)*V(:,1)' + 2.45*sqrt(D(2,2))*sin(theta)*V(:,2)' + ones(size(theta))*mean_d
            % plot(d(:,1),d(:,2),'k.');  %scatter plot with x=column 1 of d, y=column 2
            % hold on;
            % plot(semimaj(:,1),semimaj(:,2),'r','LineWidth',2);
            % plot(semimin(:,1),semimin(:,2) ,'r','LineWidth',2);
            % plot(ellipse(:,1),ellipse(:,2) ,'g','LineWidth',2);
            %
            %
            
            %% Area95_perSec     Area of the 95% confidence ellipse  per second
            % see above for detail in calcuations
            divSec = 1:sampRate:length(COPml)-1; % divide into seconds
            for i = 1:length(divSec)
                if i == length(divSec)
                    d = [COPml(divSec(i):end), COPap(divSec(i):end)];
                else
                    d = [COPml(divSec(i):divSec(i)+sampRate-1), COPap(divSec(i):divSec(i)+sampRate-1)];
                end
                cov_mtx=cov(d);
                [V,D]=eig(cov_mtx);
                c = 5.991; % for 95%. For 90%, use c = 4.605
                Area95_perSec(i) = c*pi*sqrt(D(1,1)*D(2,2));
            end
            metrics.Area95_perSec = mean(Area95_perSec);
            
            %% Area95_perSec_std Area of the 95% confidence ellipse  per second, standard deviation
            metrics.Area95_perSec_std = std(Area95_perSec);
            
            %% Hull              Hull of COP (NOT IMPLEMENTED YET)
            %     %maxr is the maximum distance in the interval; maxd the corresponding angle
            %     A = sum(0.5*maxr(a)*maxr(a+1)*sind(maxd(a+1)-maxd(a)));
            %     metrics.Hull = sum(A);
            
            %% Circ95            95% prediction circumference area
            metrics.Circ95 = pi*(mean(rad) + 1.96*std(rad))^2;
            
            %% Circ95_perSec     95% prediction circumference area per second
            divSec = 1:sampRate:length(COPml)-1; % divide into seconds
            for i = 1:length(divSec)
                if i == length(divSec)
                    r = rad(divSec(i):end);
                else
                    r = rad(divSec(i):divSec(i)+sampRate-1);
                end
                Circ95_perSec(i) = pi*(mean(r) + 1.96*std(r))^2;
            end
            metrics.Circ95_perSec = mean(Circ95_perSec);
            
            %% Circ95_perSec_std 95% prediction circumference area per second, standard deviation
            metrics.Circ95_perSec_std = std(Circ95_perSec);
            
            %% TURNi             Turns index
            % A measure used to quantify the twisting and turning of the COP
            % trajectory is the 'turn index' which is the sway- path length
            % of the normalized posturogram.
            metrics.TURNi = sum(sqrt(diff(COPap/std(COPap)).^2 + diff(COPml/std(COPml)).^2));
            
            %% Beta              mean angle deviance from AP
            beta = 90-abs(atand(diff(COPap)./diff(COPml)));
            metrics.Beta = mean(beta);
            
            %% Beta_std          mean angle deviance from AP
            metrics.Beta_std = std(beta);
            
            %% Rml               Range of COP ml displacement
            metrics.Rml = range(COPml);
            
            %% Rap               Range of COP ap displacement
            metrics.Rap = range(COPap);
            
            %% SPml              Sway Path ml
            % total displacement of sway, i.e. the 'size' or length of the
            % COP trajectory,
            metrics.SPml = sum(abs(COPml));
            % for units, might divide by the sample rate to get comparable interpretation.
            %For example:
            %         ml = 0.001*ones(30000,1); % 30 seconds, at 1000 Hz
            %         SPml = sum(abs(ml)) % = 30.0 (units = ??)
            %
            %         ml = 0.001*ones(3000,1); % 30 seconds, at 100 Hz
            %         SPml = sum(abs(ml)) % = 3.0 (units = ??)
            %
            %         ml = 0.001*ones(30000,1); % 30 seconds, at 1000 Hz
            %         SPml = sum(abs(ml))/1000 % = 0.030 (units = ?? but looks like meters)
            %
            %         ml = 0.001*ones(3000,1); % 30 seconds, at 100 Hz
            %         SPml = sum(abs(ml))/100 % = 0.030 (units = ?? but looks like meters)
            
            %% SPap              Sway Path ap
            % total displacement of sway, i.e. the 'size' or length of the
            % COP trajectory,
            metrics.SPap = sum(abs(COPap));
            % for units, might divide by the sample rate to get comparable interpretation.
            %For example:
            %     ap = 0.002*ones(30000,1); % 30 seconds, at 1000 Hz
            %     SPap = sum(abs(ap)) % = 60.0 (units = ??)
            %
            %     ap = 0.002*ones(3000,1); % 30 seconds, at 100 Hz
            %     SPap = sum(abs(ap)) % = 6.0 (units = ??)
            %
            %     ap = 0.002*ones(30000,1); % 30 seconds, at 1000 Hz
            %     SPap = sum(abs(ap))/1000 % = 0.060 (units = ?? but looks like meters)
            %
            %     ap = 0.002*ones(3000,1); % 30 seconds, at 100 Hz
            %     SPap = sum(abs(ap))/100 % = 0.060 (units = ?? but looks like meters)
            
            %% SPr               Resultant Sway Path
            % note: this is very similar to the average radial displacement (RDavg)
            % total displacement of sway, i.e. the 'size' or length of the
            % COP trajectory,
            metrics.SPr = sum(sqrt(COPap.^2+COPml.^2));
            
            %% PL                Pathlength
            % note: proportional to Vavg
            metrics.PL = sum(sqrt(diff(COPap).^2 + diff(COPml).^2));
            
            %% PLml              Pathlength, ML
            % note: proportional to Vml
            metrics.PLml = sum(abs(diff(COPml)));
            
            %% PLap              Pathlength, AP
            % note: proportional to Vap
            metrics.PLap = sum(abs(diff(COPap)));
            
            %% PLn               Normalized Pathlength
            % note: proportional to Vavg
            timeDuration = length(COPml)*1/sampRate;
            metrics.PLn = metrics.PL/timeDuration;
            
            %% PLmln             Normalized Pathlength, ML
            % note: proportional to Vml
            metrics.PLmln = metrics.PLml/timeDuration;
            
            %% PLapn             Normalized Pathlength, AP
            % note: proportional to Vap
            metrics.PLapn = metrics.PLap/timeDuration;
            
            %% Power spectral density ml
            %       Peak (Fpeakml)
            %       Mean (Fmeanml)
            %       Median (F50ml)
            %       Frequency band that contains up to 80% of the spectrum (F80ml)
            
            %
            % original Duarte Code is outdated (uses psd function):
            % nfft = round(length(COPml)/2);
            % [p,f] = psd(detrend(COPml),nfft,sampRate,nfft,round(nfft/2));
            
            % adapted from code in Schubert 2012 (J Bio Sci & Eng)
            nfft = round(length(COPml)/2);
            [p,f]=pwelch(detrend(COPml),[],[],nfft,sampRate);
            
            % from Duarte
            [m, peak] = max(p);
            area = cumtrapz(f,p);
            find50 = find(area>=.50*area(end));
            find80 = find(area>=.80*area(end));
            metrics.Fmeanml = trapz(f,f.*p)/trapz(f,p);
            metrics.Fpeakml = f(peak);
            metrics.F50ml = f(find50(1));
            metrics.F80ml = f(find80(1));
            % plot(f,p)
            
            %% Power spectral density ap
            %       Peak (Fpeakml)
            %       Mean (Fmeanml)
            %       Median (F50ml)
            %       Frequency band that contains up to 80% of the spectrum (F80ml)
            
            nfft = round(length(COPap)/2);
            [p,f]=pwelch(detrend(COPap),[],[],nfft,sampRate);
            
            % from Duarte
            [m, peak] = max(p);
            area = cumtrapz(f,p);
            find50 = find(area>=.50*area(end));
            find80 = find(area>=.80*area(end));
            metrics.Fmeanap = trapz(f,f.*p)/trapz(f,p);
            metrics.Fpeakap = f(peak);
            metrics.F50ap = f(find50(1));
            metrics.F80ap = f(find80(1));
            % plot(f,p)
            
            %% Fast Fourier Transform (FFT)
            % following matlab example: https://www.mathworks.com/help/matlab/ref/fft.html
            L = length(COPml); % length of signal
            fs = sampRate;
            
            f = fs*(0:(L/2))/L;  % Define the frequency domain (single side)
            
            % COPml
            Y = fft(detrend(COPml)); % Discrete fourier transform (using a Fast Fourier Transform algorithm)
            
            % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
            P2 = abs(Y/L); % note: normalized Fourier transform
            P1 = P2(1:L/2+1); % can give warning if not an even length
            P1(2:end-1) = 2*P1(2:end-1); % P1: the single-sided amplitude spectrum
            %ph2 = angle(Y/L)*180/pi;
            %ph1 = ph2(1:L/2+1); % phase spectrum
            FFTml = P1; % amplitudes
            
            % COPap
            Y = fft(detrend(COPap)); % Discrete fourier transform (using a Fast Fourier Transform algorithm)
            
            % Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
            P2 = abs(Y/L); % note: normalized Fourier transform
            P1 = P2(1:L/2+1); % can give warning if not an even length
            P1(2:end-1) = 2*P1(2:end-1); % P1: the single-sided amplitude spectrum
            %ph2 = angle(Y/L)*180/pi;
            %ph1 = ph2(1:L/2+1); % phase spectrum
            FFTap = P1; % amplitudes
            
            % save output FFT data
            metrics.FFTfrequencies = f;
            metrics.FFTml = FFTml;
            metrics.FFTap = FFTap;
            
            
            
            
        end
        function COPavg = averageMetrics(COP,groups, groupLabels)
            % compute average metrics based on groupings specified
            % INPUT: 
            %      COP: structure, where metrics have already been calculated 
            %      groups: Optional. cell array specifying which trials belong to which group
            %              cell array [groupNum x 1], each element containing vector of trial numbers. 
            %      groupLabels: optional. cell array [groupNum x 1] names of the groups
            % 
            % for example: Conducted 6 posturography trials, 3 eyes closed,
            % 3 eyes open. we could specify these 2 groups of trials as:
            %       groups = {[1 2 3]; [4 5 6]}; 
            %       groupLabels = {'eyes closed','eyes open'};
            % 
            
            % default values
            if nargin == 1
                groups = {[1:length(COP)]};
                groupLabels = {'All trials'};
            end
            if nargin == 2
                for i = 1:length(groups)
                    groupLabels{i} = char(i+'A'-1); % if no labels, use letters
                end
            end
            
            assert(length(groups)==length(groupLabels)); 
            
            for g = 1:length(groups)
                for i = 1:length(groups{g}) % pull trials in that group
                    metrics_trial(i) = COP(groups{g}(i)).metrics;
                end
                if length(groups{g})==1 % dont need to average if group of 1 member
                    COPavg(g).filename = groupLabels{g};
                    COPavg(g).nTrials = 1;
                    COPavg(g).metrics = metrics_trial;
                    continue;
                end
                % set up new structure:
                fnames = fieldnames(metrics_trial);
                c = cell(length(fnames),1);
                metrics=cell2struct(c,fnames);
                % fill struct with average values
                for i = 1:length(fnames)
                    if ~startsWith(fnames{i},'FFT') 
                        % use normal average scheme
                        metrics.(fnames{i}) = mean([metrics_trial.(fnames{i})]);
                    end
                end
                % use different averaging scheme for the FFT metrics
                for j = 1:length(groups{g}) % ensure FFTfrequencies vectors are the same
                    try
                    FFTfrequencies(:,j) = metrics_trial(j).FFTfrequencies;
                    catch
                        error(['FFTfrequencies in group: ' num2str(g) ' are different. May have to interpolate when averaging, but this is not implemented yet. Or the lengths of the FFT vectors are different, need to update code to accomodate these situations.'])
                    end
                end
                for j = 1:length(groups{g}) % pull FFTml and FFTap values
                    FFTml(:,j) = metrics_trial(j).FFTml;
                    FFTap(:,j) = metrics_trial(j).FFTap;
                end
                % set average FFT output
                metrics.FFTfrequencies = FFTfrequencies(:,1);
                metrics.FFTml = mean(FFTml,2);
                metrics.FFTap = mean(FFTap,2);
                
                
                % save to output
                COPavg(g).filename = groupLabels{g};
                COPavg(g).nTrials = length(groups{g});
                COPavg(g).metrics = metrics;
            end
        end
        function fig = plotMetrics(COP,order,errorType)
            % plots the posturography metrics.
             % INPUT: 
            %      COP: structure, such that COP.ml,COP.ap, COP.time exist
            %      order: Optional. array specifying how to group metrics
            %      errorType: 1 = standard error error bars, 2 = standard deviation error bars
            % for example: if order = [1 2 3; 4 5 6]; then 2 groups of 3
            % reps each
            if nargin == 1
                order = 1:length(COP);
            end
            
            if nargin == 2
                errorType = 1; % standard error is the default error bars
            end
            
            % plot options
            lineWidth = 2;
            fontSize = 12;
            
            [ngroups,nreps] = size(order);
            
            % pull reps of a group
            groupVals = cell(ngroups,1);
            for i = 1:ngroups
                groupVals{i} = [COP(order(i,:)).metrics]; % pull all reps of a group
            end
            
            % extract metric names
            metricNames = fieldnames(COP(1).metrics);
            metricNames = metricNames(1:end-3);
            
            % cycle through all metrics
            fig = figure(); fig.Name = 'posturographyMetrics';
            tabgp = uitabgroup();
            intervals = [1:6:length(metricNames)];
            tabCount = 0;
            for fn = 1:length(metricNames)
               
                % find mean and std for each group
               for g = 1:length(groupVals)
                   metric_mean(g)    = mean([groupVals{g}.(metricNames{fn})]); % mean
                   metric_std(g)     =  std([groupVals{g}.(metricNames{fn})]); % standard deviation
                   metric_stderror(g)=  std([groupVals{g}.(metricNames{fn})])/sqrt(length([groupVals{g}.(metricNames{fn})])); % standard error
               end
               
               % plot
               if any(fn==intervals)
                   tabCount = tabCount+1;
                   if fn~=intervals(end)
                       tabName = ['Metrics ' num2str(intervals(tabCount)) ' - ' num2str(intervals(tabCount+1)-1)];
                   else
                       tabName = ['Metrics ' num2str(intervals(tabCount)) ' - ' num2str(length(metricNames))];
                   end
                   tab(tabCount) = uitab(tabgp,'title',tabName,'BackgroundColor',[1 1 1]);
                   if errorType == 1; errorBarNotice = 'Error bars: standard error'; elseif errorType == 2; errorBarNotice = 'Error bars: standard deviation'; end;
                   annotation(tab(tabCount),'textbox', [0.5, 0, 0.05, 0.05], 'string', errorBarNotice,'FitBoxToText','on','HorizontalAlignment','center');
                   sp = 1; % starting subplot value
               end
               subplot(2,3,sp,'parent',tab(tabCount)); sp = sp+1;
               x = [1:length(groupVals)];
               bar(x,metric_mean)
               hold on
               if errorType == 1
                   errorbar(x,metric_mean,metric_stderror,'k.','LineWidth',lineWidth)
               elseif errorType == 2
                   errorbar(x,metric_mean,metric_std,'k.','LineWidth',lineWidth)
               else
                   error('Unknown parameter: errorType');
               end
               hold off
               title(metricNames{fn},'FontSize',fontSize,'Interpreter','none');
               xlabel('Group','FontSize',fontSize)
            end
            
        end
        function fig = plotSpectralAnalysis(COP,order,xlimits)
            % plot spectral analysis (FFT) values.
            % INPUT: 
            %      COP: structure, where metrics have already been calculated 
            %      order: Optional, depending on plot type (see below)
            %      xlimits: optional. replaces default x limits for plots
            % 
            % plot types:
            % 1: single COP, 1 figure with ML and AP side by side
            % 2: multiple COPs, 2 figures (ML,AP) with order of plots specified
            %       order: Optional. array specifying how to group metrics
            %           for example: if order = [1 2 3; 4 5 6]; then plots as
            %           6 subplots, with 2 rows, 3 columns
            % 3: multiple COPs, 2 figures (ML,AP) with order specified and
            % overlaid plots specified
            %       order: required to do overlaid plots. cell array
            %       specifying which plots will be overlaid.
            %           for example: if order = {[1 2],[3 4];[5 6],[7 8]},
            %           then plots as 4 subplots, each with 2 overlaid plots
            
            % default values
            if nargin < 3
                xlimits = [0 5]; % default range of 0-5Hz
            end
            if nargin < 2
                order = 1:length(COP);
            end
            
            % ensure FFT values exist
            for i = 1:length(COP)
                if ~all(isfield(COP(i).metrics,{'FFTfrequencies','FFTml','FFTap'}))
                    error('Must calculate all the FFT values first. Use function: posturography.calcMetrics')
                end
            end
            
            % find max y values
            ylimMaxML = 0; ylimMaxAP = 0;
            for i = 1:length(COP)
                if max(COP(i).metrics.FFTml) > ylimMaxML
                    ylimMaxML = max(COP(i).metrics.FFTml);
                end
                if max(COP(i).metrics.FFTap) > ylimMaxAP
                    ylimMaxAP = max(COP(i).metrics.FFTap);
                end
            end
            ylimMax = max([ylimMaxML ylimMaxAP]);
            
            % set plotting options
            opt.xlimits = xlimits;
            opt.ylimMax = ylimMax;
            opt.lineWidth = 2;
            opt.fontSize = 12;
            opt.unitConversion=1000;
            opt.units = 'mm';
            
            % setup figures
            fig(1) = figure(); % ML
            fig(1).PaperUnits = 'inches';
            fig(1).PaperPosition = [0 0 8.5 11];
            fig(1).Name = ['spectralAnalysis_ML'];
            fig(2) = figure(); % AP
            fig(2).PaperUnits = 'inches';
            fig(2).PaperPosition = [0 0 8.5 11];
            fig(2).Name = ['spectralAnalysis_AP'];
            
            % determine which type of plot to make
            if length(COP) == 1 % plot single COP
                delete(fig(2)); % only need one figure
                fig(1).Name = ['spectralAnalysis'];
                axML = subplot(1,2,1);
                posturography.plotSpectralAnalysis_plotter('ML',COP,opt); % outsource to plotter
                axAP = subplot(1,2,2);
                posturography.plotSpectralAnalysis_plotter('AP',COP,opt); % outsource to plotter
            else % plot multiple COPs (even overlaid)
                % configure subplot arrangements
                plot_size = size(order);
                nPlots = 1:numel(order);
                order2 = order'; % need to transpose for some reason
                
                if ~iscell(order2) % use consistent order type (overlaid vs regular)
                    order2 = num2cell(order2);
                end
                
                % plot each subplot
                for i = nPlots
                    figure(fig(1)); % ML
                    axML(i) = subplot(plot_size(1),plot_size(2),i);
                    posturography.plotSpectralAnalysis_plotter('ML',COP(order2{i}),opt); % outsource to plotter
                    
                    figure(fig(2)); % AP
                    axAP(i) = subplot(plot_size(1),plot_size(2),i);
                    posturography.plotSpectralAnalysis_plotter('AP',COP(order2{i}),opt); % outsource to plotter
                end
            end
            linkaxes([axML axAP],'xy'); % link axes so they zoom together
        end
        function saveFigure(fig,saveType,appendString)
            % inputs:
            %  fig: handle for a figure
            %  saveType: 1 for png, 2 for epsc
            %  appendString: anything you want appended to save name
            
            if nargin == 2
                appendString = '';
            elseif nargin ==3 || appendString(1) ~= '_'
                appendString = ['_' appendString];
            end
            
            % add plot directory, if it doesnt exist yet
            if ~(7==exist('plots','dir')) 
                mkdir('plots');
            end
            
            % focus on figure
            figure(fig);
            
            % set save type
            switch saveType
                case 1 % save a png
                    saveTypeCode = '-dpng';
                case 2 % save a epsc
                    saveTypeCode = '-depsc';
            end
            
            if strcmp(class(fig.Children),'matlab.ui.container.TabGroup') % check if this is a tabbed figure
                tabGroup = fig.Children;
                nChildren = length(tabGroup.Children);
                for i = 1:nChildren
                    tabGroup.SelectedTab = tabGroup.Children(i);
                    print(['plots/' fig.Name num2str(i) appendString],saveTypeCode,'-painters','-loose');
                end
            else % not a tabbed figure
                print(['plots/' fig.Name appendString],saveTypeCode,'-painters','-loose');
            end
        end
    end
    methods (Static, Access = private)
        % insert methods that can only be accessed by this class.
        function plotSpectralAnalysis_plotter(direction,COP,opt)
            % plotSpectralAnalysis outsources to this function to do the actual plotting
            
            % plotting options, must be specified.
            if nargin < 3
                % defaults
                opt.xlimits = [0 5];
                opt.ylimMax = 5;
                opt.lineWidth = 2;
                opt.fontSize = 12;
                opt.unitConversion =1000;
                opt.units = 'mm';
            end
            
            % plotting COP data in which direction?
            switch direction
                case {'ML','ml',1}
                    FFTmetric = 'FFTml';
                    FFTlabel = 'ML';
                    FFTsublabel = '_M_L';
                case {'AP','ap',2}
                    FFTmetric = 'FFTap';
                    FFTlabel = 'AP';
                    FFTsublabel = '_A_P';
            end
            
            % plot
            hold on
            for i = 1:length(COP)
                plot(COP(i).metrics.FFTfrequencies,opt.unitConversion*COP(i).metrics.(FFTmetric),'LineWidth',opt.lineWidth);
                if isfield(COP,'filename')
                    filenames{i} = COP(i).filename; % store filename
                else
                    filename{i} = char(i+'A'-1); % if no filename, use letters
                end
            end
            hold off
            
            % axes
            xlim(opt.xlimits);
            ylim([0,opt.unitConversion*opt.ylimMax]);
            
            
            % labels
            xlabel('Frequency (Hz)','FontSize',opt.fontSize)
            ylabel(['COP' FFTsublabel ' Displacement (' opt.units ')'],'FontSize',opt.fontSize)
            
            
            if length(COP) > 1 % if overlaid plots...
                legend(filenames,'Interpreter','none'); % add legend
                title(['Spectral Analysis: ' FFTlabel ' Direction'],'Interpreter','none','FontSize',opt.fontSize); % simplify title
            else
                title(['FFT, ' FFTlabel ': ' COP.filename],'Interpreter','none','FontSize',opt.fontSize); % filename as title
            end
        end
    end
end

