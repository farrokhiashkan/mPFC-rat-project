%% Wavelet Spectrogram
% Type: Morlet
% Date: '15-Jan-2020'
% Last edit: '28-May-2021'
% Ashakan Farrokhi

%=== Input:
% Signal     :  [(Trial/Channel) * Sample] 
% freq2use   :  vector contain frequency of interest
% num_cycles :  number of cycles for each frequency
% OutType    :  type of output Spectrogram ...     
% TrialMean    :  mean spectrogram accross Trials: 'Yes' or 'No'

%=== OutPut:
% Spec       :  Spectrogram of signal [Freqeuncy * Time*(Trial/Channel)]
% Comx       :  Analytic filtered signal


%% Main Function:
function [Spec,Comx,Freq2Use,num_cycles] = AFMWspectrogram(Signal, varargin)

%=== Signal Size:
SigTrial=size(Signal,1);
SigPoint=size(Signal,2);

%=== Frequency Resolution:
Fs=1000;
Freq2Use  = logspace(log10(.4),log10(50),200);
num_cycles    =  logspace(log10(3),log10(12),length(Freq2Use));

TrialMean='Yes';
PlotOption='On';

while ~isempty(varargin)
    switch lower(varargin{1})
        
        case 'fs'
            Fs = varargin{2};
            
        case 'freq2use'
            Freq2Use = varargin{2};
            
        case 'num_cycles'
            num_cycles = varargin{2};
            
        case 'trialmean'
            TrialMean = varargin{2};
            
        case 'plot'
            PlotOption = varargin{2};
            
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    varargin(1:2) = [];
end



%=== Spectrogram
if strcmpi(TrialMean,'Yes')
    Spec=zeros(length(Freq2Use),SigPoint);
else
    Spec=zeros(length(Freq2Use),SigPoint,SigTrial);
end
Comx=zeros(length(Freq2Use),SigPoint,SigTrial);

%=== wavelet and FFT parameters
ts=Fs*(num_cycles(1)/(2*pi*Freq2Use(1)));
time          = -round(4*ts)/Fs:1/Fs:round(4*ts)/Fs;
half_wavelet  = fix((length(time)-1)/2);
n_wavelet     = length(time);
n_data        = SigPoint*SigTrial;
n_convolution = n_wavelet+n_data-1;

for f=1:length(Freq2Use)
    
    % create wavelet and take FFT
    s = num_cycles(f)/(2*pi*Freq2Use(f));
    waveletSig=exp(2*1i*pi*Freq2Use(f).*time) .* exp(-time.^2./(2*(s^2)));
    
    waveletSig=sqrt(2).*(waveletSig./sum(abs(waveletSig).^2));
    wavelet_fft = fft( waveletSig', n_convolution);
    
    x=Signal';
    x=x(:);
    x = fft(x,n_convolution);
    
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*x,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    sig1 = reshape(convolution_result_fft,SigPoint,SigTrial);
    
    Comx(f,:,:)=sig1;
    sig1=abs(sig1).^2;
    
    if strcmpi(TrialMean,'Yes')==1
        Spec(f,:)=nanmean(sig1,2);
    else
        Spec(f,:,:)=sig1;
    end 
    
end %f
 

if strcmpi(PlotOption,'On')
    
    %%%%%%%%%%%%%%%%%%%
    figure
    imagesc(mean(Spec,3))
    title('Morlet Wavelet')
    ax=gca();
    ax.YDir='normal';
    ax.XDir='normal';
    Ticy=yticks;
    yticks((Ticy))
    yticklabels((Freq2Use(Ticy)))
    xlabel('Time')
    ylabel('Frequency(Hz)')
    c = colorbar;
    c.Label.String = 'instantaneous power';
    
    %%%%%%%%%%%%%%%%%
    figure
    imagesc(mean(abs(Comx),3))
    title('Morlet Wavelet')
    ax=gca();
    ax.YDir='normal';
    ax.XDir='normal';
    Ticy=yticks;
    yticks((Ticy))
    yticklabels((Freq2Use(Ticy)))
    xlabel('Time')
    ylabel('Frequency(Hz)')
    c = colorbar;
    c.Label.String = 'magnitude';
end

end
