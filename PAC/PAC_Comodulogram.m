% PAC Calculation 
% Comodulogram
% Connectivity calculation
% For T-maze data
%-----------------------------------------------------------------
% Date: '20-Apr-2020'
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
Freq_Low = [3 1 12];                                                        % [F_Phase_low StepSize F_Phase_high]
Freq_High = [15 5 100];                                                     % [F_Amp_low   StepSize F_Amp_high]
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
[~,mvl_Normal] = AF_PAC_mvlNormal(FilterDesign,Xamp,Xphs,'CatTrial','Gaurd',Gaurd,'Nbin',Nbin);
PAC = abs(mean(mvl_Normal,3));


disp('Calculate PAC_shuffled')
N_Shuffle =1;                                                             % number of shuffled data(trial shuffling)
PAC_Shuffle = zeros(size(PAC,1),size(PAC,2),N_Shuffle);
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
  
    [~,mvl_Normal] = AF_PAC_mvlNormal(FilterDesign,Xamp,Xphs,'CatTrial','Gaurd',Gaurd);
    pac = abs(mean(mvl_Normal,3));
    
    PAC_Shuffle(:,:,Sh) = pac ;  

end %Sh

%% Plot: PAC comodulogram

Center_Low = FilterDesign.Fc_Low ;
Center_High = FilterDesign.Fc_High ;

Input = (PAC -  mean(PAC_Shuffle,3) - (2*std(PAC_Shuffle,[],3)) ) ;
Input(Input<0)=0;

Input_Normalized = interp2(Input,4);
Input = interp2(PAC,4);
Input_Shuffled = interp2(mean(PAC_Shuffle,3),4);


figure(1)
imagesc(Center_High,Center_Low,Input);
colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('frequency(Hz)')
ylabel('frequency(Hz)')
colorbar
title('PAC')


figure(2)
imagesc(Center_High,Center_Low,Input_Shuffled);
colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('frequency(Hz)')
ylabel('frequency(Hz)')
colorbar
title('PAC Shuffled')


figure(3)
imagesc(Center_High,Center_Low,Input_Normalized);
colormap jet
ax=gca ;
ax.YDir='normal';

xlabel('frequency(Hz)')
ylabel('frequency(Hz)')
colorbar
title('Normalized PAC')





