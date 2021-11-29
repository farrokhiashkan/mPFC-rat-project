%% This example explained how use common frequency pattern (CFP) for classification
%% Bootstrap sampling, using replicates, was performed to evaluate the...
% performance of comparisons.
%% sensitivity analysis was proposed to evaluate the contribution of each ...
% frequency component to discrimination
%% Trial shuffling approch performed for calculating chance accuracy 
%%
clear
clc
close all hidden
load ('SampleData.mat')

%===parameter for CFP 
% alphaSet=[.01, .1, .2, .4, .5, 1, 2, 5, 10, 100];
alphaSet=.1;
m=3;

%% %%%%%%%%%%%%%%Classification%%%%%%%%%%%%%%%

Train_Percent = .8;                                                         % percent of data use for training
N_Rand_Sampleing= 10;                                                       % number of replicates for bootstrap
NTr1 = round(Train_Percent * size(Class1,1));
NTr2 = round(Train_Percent * size(Class2,1));
Indices1 = 1:size(Class1,1);
Indices2 = 1:size(Class2,1);

mgs='-';
fprintf(['   Accuracy=',mgs]);


HeatMap_Total = zeros(1,size(Class1,3));
PERF = nan(N_Rand_Sampleing,1);
for K=1:N_Rand_Sampleing
    
    Tr1_Index =  randperm(size(Class1,1),NTr1) ;
    Tr2_Index =  randperm(size(Class2,1),NTr2) ;
    
    Te1_Index =  Indices1 ; Te1_Index(Tr1_Index)=[];
    Te2_Index =  Indices2 ; Te2_Index(Tr2_Index)=[];
    
    
    %=== Train
    Train1 = Class1(Tr1_Index,:,:);
    Train2 = Class2(Tr2_Index,:,:);
    
    %=== Test
    Test1 = Class1(Te1_Index,:,:);
    Test2 = Class2(Te2_Index,:,:);
    
    %=== CSP
    [W,alpha,C1,C2]=AF_CFP(Train1,Train2,m,alphaSet);
    
    %=== Training 1
    Ftrain1=nan(size(Train1,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Train1,1)
        x=squeeze(Train1(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftrain1(Trial,:)=f;
    end %Trial
    
    %=== Training 2
    Ftrain2=nan(size(Train2,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Train2,1)
        x=squeeze(Train2(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftrain2(Trial,:)=f;
    end %Trial
    
    %=== Test 1
    Ftest1=nan(size(Test1,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Test1,1)
        x=squeeze(Test1(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftest1(Trial,:)=f;
    end %Trial
    
    %=== Test 2
    Ftest2=nan(size(Test2,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Test2,1)
        x=squeeze(Test2(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftest2(Trial,:)=f;
    end %Trial
    clear x y f Trial;
    
    Ftrain = [Ftrain1;Ftrain2];
    L_Train = [zeros(size(Ftrain1,1),1); ones(size(Ftrain2,1),1)];
    
    Ftest = [Ftest1;Ftest2];
    L_Test = [zeros(size(Ftest1,1),1); ones(size(Ftest2,1),1)];
    
    SVMStruct = fitcsvm(Ftrain,L_Train);
    [pred,~] = predict(SVMStruct,Ftest);
    PERF(K)= mean(pred==L_Test);
    
    %=== Heat map:
    HeatMap=nan(1,size(W,1));
    for fi = 1:size(HeatMap,2)
        TestHeat = zeros(1,size(HeatMap,2));
        TestHeat =  TestHeat * (W.^2);
        [~,Score1] = predict(SVMStruct,TestHeat);
        
        TestHeat = zeros(1,size(HeatMap,2));
        TestHeat(1,fi) = TestHeat(1,fi)+ 10^(-1) ;
        TestHeat =  TestHeat * (W.^2);
        [~,Score2] = predict(SVMStruct,TestHeat);
        HeatMap(1,fi) = abs( Score2(1)-Score1(1) ) /(10^(-2));
    end
    HeatMap_Total= HeatMap_Total + (HeatMap-min(HeatMap,[],'all')) / (max(HeatMap,[],'all')-min(HeatMap,[],'all'));
    
    fprintf(repmat('\b',1,length(mgs)))
    mgs = sprintf([num2str(fix(nanmean(PERF)*10000)/100), ' (',num2str(K),'/',num2str(N_Rand_Sampleing),')','  alpha=',num2str(alpha)]);
    fprintf(mgs)
end %K
HeatMap_Total = HeatMap_Total/K;
HeatMap_Total = (HeatMap_Total-min(HeatMap_Total,[],'all')) / (max(HeatMap_Total,[],'all')-min(HeatMap_Total,[],'all'));

PerfTOT = mean(PERF);
PerfTOT_Seg{1}.Perf = PerfTOT;

W_Total{1} = HeatMap_Total;

fprintf(repmat('\b',1,length(mgs)))
fprintf(num2str(fix(PerfTOT*10000)/100))

figure(1)
area(HeatMap_Total,'FaceColor','b')
xlabel('Frequency(Hz)')
ylabel('dScore/dPower')
sgtitle('BT2','Interpreter', 'none')
drawnow


%% Chance Accuracy
N_Shuffle = N_Rand_Sampleing;                                               % number of shuffled dataset

Class_Cat = cat(1,Class1,Class2);
Ind = randperm(size(Class_Cat,1),size(Class_Cat,1));
Class_Cat = Class_Cat(Ind,:,:);

mgs='-';
fprintf(['   Chance_Accuracy= ',mgs]);

Perf_Shuffle = nan(N_Shuffle,1);
L_Class = size(Class1,1);
for Sh= 1:N_Shuffle
    
    Ind = randperm(size(Class_Cat,1),size(Class_Cat,1));
    Class_Cat_Shuffle = Class_Cat(Ind,:,:);
    Class1_Shuffle = Class_Cat_Shuffle(1:L_Class,:,:) ;
    Class2_Shuffle = Class_Cat_Shuffle(end-L_Class+1 : end,:,:) ;
    
    Tr1_Index =  randperm(size(Class1_Shuffle,1),NTr1) ;
    Tr2_Index =  randperm(size(Class2_Shuffle,1),NTr2) ;
    
    Te1_Index =  Indices1 ; Te1_Index(Tr1_Index)=[];
    Te2_Index =  Indices2 ; Te2_Index(Tr2_Index)=[];
    
    %=== Train
    Train1 = Class1_Shuffle(Tr1_Index,:,:);
    Train2 = Class2_Shuffle(Tr2_Index,:,:);
    
    %=== Test
    Test1 = Class1_Shuffle(Te1_Index,:,:);
    Test2 = Class2_Shuffle(Te2_Index,:,:);
    
    
    [W,~,~,~]=AF_CFP(Train1,Train2,m,alphaSet);
    
    %=== Training 1
    Ftrain1=nan(size(Train1,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Train1,1)
        x=squeeze(Train1(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftrain1(Trial,:)=f;
    end %Trial
    
    %=== Training 2
    Ftrain2=nan(size(Train2,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Train2,1)
        x=squeeze(Train2(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftrain2(Trial,:)=f;
    end %Trial
    
    %=== Test 1
    Ftest1=nan(size(Test1,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Test1,1)
        x=squeeze(Test1(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftest1(Trial,:)=f;
    end %Trial
    
    %=== Test 2
    Ftest2=nan(size(Test2,1),size(W,2));
    clear x; clear y; clear f;
    for Trial = 1:size(Test2,1)
        x=squeeze(Test2(Trial,:,:));
        x=x*W;
        f=var(x);
        Ftest2(Trial,:)=f;
    end %Trial
    clear x y f Trial;
    
    Ftrain = [Ftrain1;Ftrain2];
    L_Train = [zeros(size(Ftrain1,1),1); ones(size(Ftrain2,1),1)];
    
    Ftest = [Ftest1;Ftest2];
    L_Test = [zeros(size(Ftest1,1),1); ones(size(Ftest2,1),1)];
    
    SVMStruct = fitcsvm(Ftrain,L_Train);
    [pred,~] = predict(SVMStruct,Ftest);
    Perf_Shuffle(Sh)= mean(pred==L_Test);
    
    fprintf(repmat('\b',1,length(mgs)))
    mgs = sprintf([num2str(fix(nanmean(Perf_Shuffle(1:Sh))*10000)/100), ' (',num2str(Sh),'/',num2str(N_Shuffle),')']);
    fprintf(mgs)
    
end %Sh


%=== Statistical Test:
P_Value(1) = ranksum(PERF,Perf_Shuffle);
fprintf(['   P_Value=',num2str(P_Value(1)),'\n'])

