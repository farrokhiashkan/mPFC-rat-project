% Last edit: '30-Nov-2021'
% Ashakan Farrokhi
% IUST

% Filter design:
% find filters coefficient used for PAC calculation 
% Low frequency filters and High frequency filters
% approach : Fix >> fix bandwidth for Low and High frequency bands
% approach : Addaptive >> fix bandwidth for Low and addaptive bands for High frequency components
% for 'Fix' approach High frequency components filtered using 5 as bandwidth

%Input:
% Fs: sampling frequency
% Mode: approach >> 'Fix' or 'Addaptive'

%Output:
% FilterDesign: contain filters parameter

% Example:
% Freq_Low = [3 1 20];
% Freq_High = [30 5 200];
% Mode = 'Addaptive';
% Fs=1000;
% FilterDesign = AF_FilterBank2(Fs,Mode,Freq_Low,Freq_High);

%Note: for Freq_Low : if second element set 'nan' >> one filter defined between first and third elements...
%>> for example : [3 nan 10] >> only one filter desined between 3 and 10 (Hz) 

%%
function FilterDesign=AF_FilterBank(Fs,Mode,Freq_Low,Freq_High)

Order=3;
%Low frequency component:
% flmin=3; %Hz
% flmax=20; %Hz
flmin = Freq_Low(1);
BLStep = Freq_Low(2);
flmax = Freq_Low(3);

if isnan(BLStep)
    Bands_Low = [flmin flmax];
    Fc_Low = mean(Bands_Low);
else
    Fc_Low=flmin:BLStep:flmax ;
    Bands_Low=[Fc_Low-1; Fc_Low+1 ]';
end

%High frequency component:
% BHStep=5; %Hz
% Fc_High(:,1) = 30:BHStep:200;
fhmin = Freq_High(1);
BHStep = Freq_High(2);
fhmax = Freq_High(3);
Fc_High(:,1) = fhmin:BHStep:fhmax;


if strcmpi(Mode,'Fix')==1
    Band_High_baseLow=5*ones(size(Bands_Low,1),1); %Hz
elseif strcmpi(Mode,'Addaptive')==1
    Band_High_baseLow= 1*(Bands_Low(:,2));
end

%design filter:
Coeff_b_Low=zeros(size(Bands_Low,1),2*Order+1);
Coeff_a_Low=zeros(size(Bands_Low,1),2*Order+1);

Coeff_b_High=cell(size(Bands_Low,1),length(Fc_High));
Coeff_a_High=cell(size(Bands_Low,1),length(Fc_High));
Bands_High=cell(size(Bands_Low,1),1);
for Lband=1:size(Bands_Low,1)
    [Coeff_b_Low(Lband,:), Coeff_a_Low(Lband,:)]=butter(Order,[Bands_Low(Lband,1) Bands_Low(Lband,2)]/(Fs/2));
    
    Bands_High{Lband,1}=[Fc_High-Band_High_baseLow(Lband),Fc_High+Band_High_baseLow(Lband)];
    for Hband=1:length(Fc_High)
        [Coeff_b_High{Lband,Hband}, Coeff_a_High{Lband,Hband}]=butter(Order,[Bands_High{Lband,1}(Hband,1) Bands_High{Lband,1}(Hband,2)]/(Fs/2));
    end %Hband
    
end %Lband


FilterDesign.Coeff_b_High=Coeff_b_High;
FilterDesign.Coeff_a_High=Coeff_a_High;
FilterDesign.Fc_High=Fc_High;
FilterDesign.Coeff_b_Low=Coeff_b_Low;
FilterDesign.Coeff_a_Low=Coeff_a_Low;
FilterDesign.Fc_Low=Fc_Low;


end