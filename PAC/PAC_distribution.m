% PAC Calculation (dsitribution)

% calculate the distribution of high freqeuncy signal components on...
% different phases of a specific lower frequency band.

% Connectivity calculation
% For T-maze data
%-----------------------------------------------------------------
% Last edit: '29-Nov-2021'
% Ashkan Farrokhi
% IUST
%% Initialization:
clear;
clc;
close all;
set(0,'DefaultFigureWindowStyle','normal')  %  'normal' or 'docked'
%% Load data
% load simulated LFP (with PAC)
% this data simulated using "AFGenerateLFPpac" function
load LFP_Simulated_Signal
Signal = Report.Signal.Xnew;                                                
%% Analysis:

%=== Design Filter Bank;
Freq_Low = [3 nan 12];                                                      % [F_Phase_low StepSize F_Phase_high]
Freq_High = [15 5 120];                                                     % [F_Amp_low   StepSize F_Amp_high]
Fs=1000;                                                                    % Frequency sampling (Hz)
Mode='Addaptive';
FilterDesign=AF_FilterBank(Fs,Mode,Freq_Low,Freq_High);

%=== Analysis:
clc
disp('Calculate PAC')

Xamp = Signal;                                                              %signal used for phase component
Xphs = Signal;                                                              %signal used for Amplitude component
Gaurd=100;                                                                  %Gaurd(not included in PAC). Signal:[t1-Gaurd t2+Gaurd]
Nbin=18;                                                                    %when set the Histogram method also performed 
[PAC_Hist,~] = AF_PAC_mvlNormal(FilterDesign,Xamp,Xphs,'CatTrial','Gaurd',Gaurd,'Nbin',Nbin);
Hist = squeeze(PAC_Hist);


disp('Calculate PAC_shuffled')
N_Shuffle =100;                                                             % number of shuffled data(trial shuffling)
Hist_Shuffle = zeros(size(Hist,1),size(Hist,2),N_Shuffle);
for Sh = 1:N_Shuffle
    
    clc;
    
    disp('Calculate PAC : completed')
    disp('Calculate PAC_shuffled')
    
    Mgs = sprintf(['N_Shuffle =',num2str(Sh),'/',num2str(N_Shuffle),'   ']);
    fprintf(Mgs)
    
    
    Ind_Trial_Phs = randperm(size(Signal,2)) ;
    Ind_Trial_Amp = randperm(size(Signal,2)) ;
    
    Xamp = Signal(:,Ind_Trial_Amp) ;
    Xphs = Signal(:,Ind_Trial_Phs) ;
    
    [hist,~] = AF_PAC_mvlNormal(FilterDesign,Xamp,Xphs,'CatTrial','Gaurd',Gaurd,'Nbin',Nbin);
    pac = squeeze(hist);
    
    Hist_Shuffle(:,:,Sh) = pac ;
    
end %Sh


%% Plot

close all;
Center_Low = linspace(0,360,18+1);
Center_Low=Center_Low(1:end-1)+(Center_Low(2:end)-Center_Low(1:end-1))/2;
Center_High = FilterDesign.Fc_High ;

Input = (Hist -  mean(Hist_Shuffle,3) - (2*std(Hist_Shuffle,[],3)) ) ;
Input(Input<0)=0;
Input = Input';


Input_Normalized = interp2(Input,4);
Input = interp2(Hist,4)';
Input_Shuffled = interp2(mean(Hist_Shuffle,3),4)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
imagesc(Center_Low,Center_High(1:end),Input);

colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('Theta Phase(degree)')
ylabel('frequency(Hz)')
colorbar
title('PAC-distribution')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2)
imagesc(Center_Low,Center_High(1:end),Input_Shuffled);

Center_Low = linspace(0,360,18+1);
Center_Low=Center_Low(1:end-1)+(Center_Low(2:end)-Center_Low(1:end-1))/2;
Center_High = FilterDesign.Fc_High ;

colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('Theta Phase(degree)')
ylabel('frequency(Hz)')
colorbar
title('PAC Shuffled')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
imagesc(Center_Low,Center_High(1:end),Input_Normalized);

Center_Low = linspace(0,360,18+1);
Center_Low=Center_Low(1:end-1)+(Center_Low(2:end)-Center_Low(1:end-1))/2;
Center_High = FilterDesign.Fc_High ;

colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('Low frequency component Phase(degree)')
ylabel('frequency(Hz)')
colorbar
title('Normalized PAC')
