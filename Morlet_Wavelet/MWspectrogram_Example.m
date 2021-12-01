%% Wavelet Spectrogram

%=== Example:
close all;
clear;
Fs=1000;

f1=60;
f2=15;
f3=6;

t=0:1/Fs:1-.001;
t2=0:1/Fs:2-.001;
Signal1(1,:)=12*sin(2*pi.*t*60);
Signal2(1,:)=6*sin(2*pi.*t*15);
Signal3(1,:)=10*sin(2*pi.*t2*6);

Signal=cat(2,Signal1,Signal2)+Signal3;

Freq2Use  = logspace(log10(1),log10(500),100); % 4-30 Hz in 15 steps
num_cycles    =  logspace(log10(3),log10(12),length(Freq2Use));
[Spec,Comx,Freq2Use,num_cycles] = AFMWspectrogram(Signal,'TrialMean','No','plot','On','fs',1000,'Freq2Use',Freq2Use,'num_cycles',num_cycles);

Ind_f1 = dsearchn(Freq2Use',f1);
Ind_f2 = dsearchn(Freq2Use',f2);
Ind_f3 = dsearchn(Freq2Use',f3);


X_f1_filt = real(Comx(Ind_f1,:));
X_f2_filt = real(Comx(Ind_f2,:));
X_f3_filt = real(Comx(Ind_f3,:));



figure
subplot(4,1,1), plot(Signal)
title('main signal')

subplot(4,1,2), plot(X_f1_filt,'b')
title('filtered signal around f1')
subplot(4,1,3), plot(X_f2_filt,'r')
title('filtered signal around f2')
subplot(4,1,4), plot(X_f3_filt,'k')
title('filtered signal around f3')




