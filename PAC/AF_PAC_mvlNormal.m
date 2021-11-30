% this function claculate the phase amplitude coupling (PAC) on input signals
% Last edit: '30-Nov-2021'
% Ashakan Farrokhi
% IUST

%=== Input:
% this function calculate PAC between related columns from Xamp and Xphs
% FilterDesign: designed filters parameter, this calculated using AF_FilterBank function
% Xamp     :  [Sample * Trial], signal for amplitude 
% Xphs   :  [Sample * Trial], signal for phase 

% varargin:
% Gaurd: gaurd removed after calculating instantaneuos phase and amp...
%for exapmle if Gaurd=100 and signal has 700 samples, then PAC calculated
%on 500 sample (101-600);

% Nbin: number of bin. if set, normalized histogram of amplitude calculate on bined phases. 
% CatTrial: trials concatenate along sample. then coupling calculate 

%=== OutPut:
% PAC_Hist     : coupling calculated using Histogram approach (dsitribution of amplitude component on bined phase component)
%[Nbin * PhaseComponent * AmpComponent * (Trial or 1)]

% mvl_Normal   : coupling calculated using normalized MVL approach(Özkurt and Schnitzler, 2011)
%[PhaseComponent * AmpComponent * (Trial or 1)]


%% Main Function:
function [PAC_Hist,mvl_Normal] = AF_PAC_mvlNormal(FilterDesign,Xamp,Xphs,varargin)

CatTrial=[];
Nbin=[];
while ~isempty(varargin)
    switch lower(varargin{1})
        
        case 'gaurd'
            Gaurd = varargin{2};
            Del_Vargin=2;
        case 'nbin'
            Nbin = varargin{2};
            Del_Vargin=2;
        case 'cattrial'
            CatTrial = 'CatTrial';
            Del_Vargin=1;
            
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
    
    if Del_Vargin==2
        varargin(1:2) = [];
    elseif Del_Vargin==1
        varargin(1) = [];
    end
end



%=== Filter Coefficients:
Coeff_b_High=FilterDesign.Coeff_b_High;
Coeff_a_High=FilterDesign.Coeff_a_High;
Fc_High=FilterDesign.Fc_High;
Coeff_b_Low=FilterDesign.Coeff_b_Low;
Coeff_a_Low=FilterDesign.Coeff_a_Low;
Fc_Low=FilterDesign.Fc_Low;

if strcmpi('CatTrial',CatTrial)
    mvl_Normal = nan(length(Fc_Low),length(Fc_High),1);
else
    mvl_Normal = nan(length(Fc_Low),length(Fc_High),size(Xamp,2));
end

%=== define Phase bin for theta(0:360 deg)
PAC_Hist = 'not calculated';
if ~isempty(Nbin)
    % Nbin=18;
    Phase_Bins=linspace(0,2*pi,Nbin+1); %define phase Bins
    if strcmpi('CatTrial',CatTrial)
        PAC_Hist = nan(Nbin,length(Fc_Low),length(Fc_High),1);
    else
        PAC_Hist = nan(Nbin,length(Fc_Low),length(Fc_High),size(Xamp,2));
    end
end %~isempty(Nbin)

for band_Low=1:length(Fc_Low)
    %=== filter Phs component
    x=filtfilt(Coeff_b_Low(band_Low,:),Coeff_a_Low(band_Low,:),Xphs);
    Phs=angle(hilbert(x));
    Phs = Phs(Gaurd+1:end-Gaurd,:);
    if strcmpi('CatTrial',CatTrial)
        Phs = Phs(:);
    end
    clear x;
    
    %=== wrapping phase between [0 2Pi]
    Phs(Phs<0) = Phs(Phs<0) + (2*pi) ;
    
    for band_High=1:length(Fc_High)
        %=== filter Amp component
        y=filtfilt(Coeff_b_High{band_Low,band_High},Coeff_a_High{band_Low,band_High},Xamp);
        Amp = abs(hilbert(y));
        Amp = Amp(Gaurd+1:end-Gaurd,:);
        if strcmpi('CatTrial',CatTrial)
            Amp = Amp(:);
        end
        clear y;
        
        %%======PAC calculation:	MVL_Normal, (Özkurt and Schnitzler, 2011)
        mvl=(Amp.*exp(1i*Phs));
        N = size(mvl,1);
        mvl = nansum(mvl) ./ sqrt( (sum(Amp.^2)) * N );
        mvl_Normal(band_Low,band_High,:) =(mvl);
        clear mvl N
        
        
        %%======PAC calculation:	Histogram(low frequency,high frequency)
        if ~isempty(Nbin)            
            Z=nan(Nbin,size(Phs,2));
            for bin=1:Nbin    %phase bin
                IND=nan(size(Phs));
                
                [ind_Row, ind_Col]=find(Phase_Bins(bin) <=Phs & Phs<= Phase_Bins(bin+1));
                ind= sub2ind(size(Phs),ind_Row,ind_Col);
                IND(ind)=1;
                
                R=IND.*Amp;
                Z(bin,:)= nanmean(R);
            end
            P = Z./nansum(Z);
            P(isnan(P))=0;
            P = permute(P,[1,3,2]);
            PAC_Hist(:,band_Low,band_High,:)=P;

            clear P            
        end %~isempty(Nbin)
             
        
    end %band_High
end %band_Low
clear Amp Phs
end