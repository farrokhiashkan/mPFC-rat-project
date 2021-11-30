% this function generate a simulated brain signal using Brownian Noise with...
% with phase amplitude coupling between user defined frequency components...
% and add noise to signal for different SNR

% first date: '26-Dec-2019'
% last edit : '30-Nov-2021'
% Ashkan Farrokhi
% IUST

% "#" >> this symbol in front of comments means that user can modify that line parameter

% Description: Simulate LFP signal: 
%>> Brain signal simulation with Brownian Noise
%>> Method of signal decomposition:
% [1]: using filter to extract components.
% [2]: add in future
%>> Method of PAC generation:
% [1]: Low frequency peaks.
% [2]: phase(instant frequency) of Low frequency.

% Report: the simulate signal saved in "LFP_Simulated_Signal.mat"
% only first trial plot during simulation

%% Initialization:
clear;
clc; 
close all;

%=== Signal spec:
Fs=1000;                                                                    %# frequency sampling(Hz) of the signal
L_Trial=2000;                                                               %# number of sample (for each trial)
N_Trial=20;                                                                 %# number of trials

%=== Brownian Noise:
cn = dsp.ColoredNoise('Color','brown','SamplesPerFrame',L_Trial,'NumChannels',N_Trial);
X= cn();
X= X-mean(X); %Zero mean 

%=== PAC spec:
Freq_Low=[5 8];      %Theta(Hz)                                             %# phase component frequency band
Freq_High=[40 50];   %Low Gamma(Hz)                                         %# Amplitude component frequency band

Gain_Amp=3;                                                                 %# Gain for phase component
Gain_Phs=3;                                                                 %# Gain for Amplitude component
SNR_add=inf;                                                                %# SNR: for adding noise
CouplingGain=.8;                                                            %# Gain of coupling (0: non coupling ,>0:stronger coupling)
%=== User Choose Method

%Method of signal decomposition:
DecMethod=1;                                                                %#  =1:Band pass Filtering or
                                                                            %#  =2: this will Add in future


%Method of PAC simulation:
PACMethod=1;                                                                %#  =1:Low Frequency Peaks or
                                                                            %#  =2:(phase)instant frequency of Low frequency

%% Signal Decomposition:
if DecMethod==1
    %===[1]: Band pass Filtering
    PlotON=0;
    [X_Phs,X_Amp,X_Oth]= FilterDec(X,Freq_Low,Freq_High,Fs,PlotON);
elseif DecMethod==2
    %===[2]: ---
end

%% Add PAC to Signal:
if PACMethod==1
    %===[1] Method: Low Frequency Peaks
    [X_Amp_Modulate,X_Phs]=ModulationPeaks(X_Phs,X_Amp,Gain_Phs,Gain_Amp,CouplingGain);
elseif PACMethod==2
    %===[2] (phase)instant frequency of Low frequency:
    [X_Amp_Modulate,X_Phs]=ModulationIP(X_Phs,X_Amp,Gain_Phs,Gain_Amp,CouplingGain);
end

Xnew = X_Amp_Modulate+X_Oth+X_Phs ;
%add white noise
Xnew_plusNoise = awgn(Xnew,SNR_add,'measured');


% Plot:
for Trial=1:N_Trial
    
    if Trial==1
        figure(Trial);     clf;
        subplot(6,1,1), plot(X(:,Trial))
        title('Raw Signal')
        title(['Trial=', num2str(Trial)])
        
        subplot(6,1,2), plot(X_Phs(:,Trial)),  hold on
        plot(angle(hilbert(X_Phs(:,Trial))),'Color','g'),  hold off
        title(strcat('Phase Signal [', num2str(Freq_Low(1)), '-', num2str(Freq_Low(2)), ' Hz]'))
        
        subplot(6,1,3), plot(X_Amp(:,Trial)),  hold on
        ax3=gca();
        %     plot(X_Amp_Modulate,'--r')
        title(strcat('Amplitude Signal [', num2str(Freq_High(1)), '-', num2str(Freq_High(2)), ' Hz]'))
        
        subplot(6,1,4), plot(X_Amp_Modulate(:,Trial))
        title('Coupled Amp Signal')
        ax4=gca();
        
        Lim3=ax3.YLim;
        Lim4=ax4.YLim;
        ax3.YLim=[min(Lim3(1),Lim4(1)) max(Lim3(2),Lim4(2))];
        ax4.YLim=[min(Lim3(1),Lim4(1)) max(Lim3(2),Lim4(2))];
    end %plot:trial
    
    Phs=angle(hilbert(X_Phs(:,Trial)));
    AmpM=abs(hilbert(X_Amp_Modulate(:,Trial)));
    
    
    Nbin=100;
    Phase_Bins=linspace(0,2*pi,Nbin+1); %define phase Bins
    Phs(Phs<0)=Phs(Phs<0)+(2*pi);
    
    
    Amp=abs(hilbert(X_Amp(:,1)));
    Z=zeros(Nbin,size(Phs,2));
    for bin=1:Nbin    %phase bin
        IND=NaN(size(Phs,1),size(Phs,2));
        
        [ind_Row, ind_Col]=find(Phase_Bins(bin) <=Phs & Phs< Phase_Bins(bin+1));
        ind= sub2ind(size(Phs),ind_Row,ind_Col);
        IND(ind)=1;
        
        R=IND.*Amp(:,:);
        Z(bin,:)= nanmean(R);
    end
    Z=Z./max(Z);
    
    if Trial==1
        subplot(6,1,5), histogram('BinEdges',rad2deg(Phase_Bins(1:end)),'BinCounts',Z)
        title('Phase amplitude plot')
        axis tight
    end %plot:trial
    
    
    Z=zeros(Nbin,size(Phs,2));
    for bin=1:Nbin    %phase bin
        IND=NaN(size(Phs,1),size(Phs,2));
        
        [ind_Row, ind_Col]=find(Phase_Bins(bin) <=Phs & Phs< Phase_Bins(bin+1));
        ind= sub2ind(size(Phs),ind_Row,ind_Col);
        IND(ind)=1;
        
        R=IND.*AmpM(:,:);
        Z(bin,:)= nanmean(R);
    end
    Z=Z./max(Z);
    
    if Trial==1
        subplot(6,1,6), histogram('BinEdges',rad2deg(Phase_Bins(1:end)),'BinCounts',Z)
        title('Phase amplitude plot')
        axis tight
    end %plot:trial
    
    
end


%%

Report.Signal.Data=X;
Report.Signal.X_Phs=X_Phs;
Report.Signal.X_Amp=X_Amp;
Report.Signal.X_Oth=X_Oth;
Report.Signal.X_Amp_Modulate=X_Amp_Modulate;
Report.Signal.Xnew=Xnew;
Report.Signal.Xnew_plusNoise=Xnew_plusNoise;


Report.Propety.Gain_Amp=Gain_Amp;
Report.Propety.Gain_Phs=Gain_Phs;
Report.Propety.SNR_add=SNR_add;
Report.Propety.CouplingGain=CouplingGain;
Report.Propety.Freq_Low=Freq_Low;
Report.Propety.Freq_High=Freq_High;
Report.Propety.Fs=Fs;


Report.Propety.DecMethod=DecMethod;
Report.Propety.PACMethod=PACMethod;
Report.Propety.Readme={'Method of signal decomposition:','=1:Band pass Filtering','=2:--';
    'Method of PAC simulation:','=1:Low Frequency Peaks','=2:(phase)instant frequency of Low frequency'};   


save('LFP_Simulated_Signal.mat','Report')
disp('signal saveds as : "LFP_Simulated_Signal.mat"')

%% Functions:

%[1]: Band pass Filtering
function [X_Phs,X_Amp,X_Oth]= FilterDec(X,Freq_Low,Freq_High,Fs,PlotON)
Order=3;

%phase component:
[b,a]= butter(Order,[Freq_Low(1) Freq_Low(2)]/(Fs/2)) ;

X_Phs= filtfilt(b,a,X);

%amplitude component:
[b,a]= butter(Order,[Freq_High(1) Freq_High(2)]/(Fs/2)) ;
X_Amp= filtfilt(b,a,X);

%Other part of signal:
X_Oth= X-(X_Amp+X_Phs);


%plot:
if PlotON==1
    figure, freqz(b,a,[],Fs)
    figure, pwelch(mean(X_Phs,2),[],[],[],Fs)
    figure , freqz(b,a,[],Fs)
    figure, pwelch(mean(X_Amp,2),[],[],[],Fs)
    figure, pwelch([mean(X,2),mean(X_Phs,2),mean(X_Amp,2),mean(X_Oth,2)],[],[],[],Fs)
    figure, plot([mean(X,2),mean(X_Phs,2),mean(X_Amp,2),mean(X_Oth,2)])
end

end %End Function: FilterDec 



%[2]: --
%...



%[1] Method: Low Frequency Peaks
function [X_Amp_Modulate,X_Phs_new]=ModulationPeaks(X_Phs,X_Amp,Gain_Phs,Gain_Amp,CouplingGain)

%Window:
L_Wind=100;  % Window length(ms)
I=CouplingGain;    % Coupling Gian [0-->inf]
WHann = I.*(hann(L_Wind+1))+1;

X_Amp_Modulate=zeros(size(X_Amp));
for Trial=1:size(X_Amp,2)
    
    %find phase component peaks:
    CouplingPhase=0;  %radian
    phase = (angle(hilbert(X_Phs(:,Trial))))-CouplingPhase;
    Locs=find(sign(phase(2:end))-sign(phase(1:end-1))~=0 );
    Locs(Locs==1)=[];
    Locs(Locs==length(phase))=[];
            
    ind=abs(phase(Locs+1)-phase(Locs-1))>1;
    Locs(ind)=[];
    
    X_Amp_Ctemp=X_Amp(:,Trial);
    for i=1:length(Locs)
        try
            X_Amp_Ctemp(Locs(i)-L_Wind/2:Locs(i)+L_Wind/2) = X_Amp_Ctemp(Locs(i)-L_Wind/2:Locs(i)+L_Wind/2).*WHann;
        catch
        end
    end
    X_Amp_Modulate(:,Trial)=Gain_Amp*X_Amp_Ctemp*(mean(abs(hilbert(X_Amp(:,Trial))))/mean(abs(hilbert(X_Amp_Ctemp))));

end %Trial

% X_Phs_new=Gain_Phs*X_Phs;
X_Phs_new=X_Phs+Gain_Phs*X_Phs;  %%%%%%%%%%%%%%%%%%%%%%%%??????????????????????????????

end %End ModulationPeaks


%[2] Method: (phase)instant frequency of Low frequency:
function [X_Amp_Modulate,X_Phs_new]=ModulationIP(X_Phs,X_Amp,Gain_Phs,Gain_Amp,CouplingGain)
C=CouplingGain;
X_Amp_Modulate=zeros(size(X_Amp));
for Trial=1:size(X_Amp,2)
    A=abs(hilbert(X_Amp(:,Trial)));
    Phs_High=angle(hilbert(X_Amp(:,Trial)));
    
    phase = (angle(hilbert(X_Phs(:,1))));
    
    a=.18;
    Pt=1-C*cos(2*pi*a*phase-(pi));
    A=Gain_Amp*(mean(A)*Pt)/mean(Pt);
    X_Amp_Modulate(:,Trial)=real(A.*exp(1i*Phs_High));
    
end %Trial
% X_Phs_new=Gain_Phs*X_Phs;
X_Phs_new=X_Phs+Gain_Phs*X_Phs;  %%%%%%%%%%%%%%%%%%%%%%%%??????????????????????????????

end %End ModulationIP

