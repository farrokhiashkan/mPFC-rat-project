% This function calculate frequency filters based of CFP algorithm
% Last edit: '29-Nov-2021'
% Ashakan Farrokhi

%=== Input:
% frequency decomposed inputs, Class1 and Class2, should feed to this function
% Class1     :  [Trial * Sample * Frequency]
% Class2   :    [Trial * Sample * Frequency]
% m :  number of pairs returned at result
% alphaSet    :  vector >> contain candidate regularization parameters (L2) , default=.01 

%=== OutPut:
% WFinal     : frequency filters
% alpha      :  parameter used for regularization
% C1         :  estimated covariance matrix for calss1
% C2         :  estimated covariance matrix for calss2


%% Main Function:
function [WFinal,alpha,C1,C2] = AF_CFP(Class1,Class2,m,alphaSet)

% determine alpha if not set
if isempty(alphaSet)
    alpha = .1;
end

if length(alphaSet)>2
    N_Fold= 5;
    
    Indices1 = crossvalind('Kfold', size(Class1,1), N_Fold);
    Indices2 = crossvalind('Kfold', size(Class2,1), N_Fold);
    PERF = nan(N_Fold,length(alphaSet));
    for K=1:N_Fold
        
        %=== Train
        Train1 = Class1(Indices1~=K,:,:);
        Train2 = Class2(Indices2~=K,:,:);
        
        %=== Test
        Test1 = Class1(Indices1==K,:,:);
        Test2 = Class2(Indices2==K,:,:);
        
        Nt1 = size(Train1,1);
        Nt2 = size(Train2,1);
        
        C1=0; C2=0;
        for Trial = 1:max(Nt1,Nt2)
            
            if Trial<=Nt1
                x=squeeze(Train1(Trial,:,:));
                x = x- mean(x,1);
                c=x'*x/trace(x'*x);
                C1=C1+c;
            end %if
            
            if Trial<=Nt2
                x=squeeze(Train2(Trial,:,:));
                x = x- mean(x,1);
                c=x'*x/trace(x'*x);
                C2=C2+c;
            end %if
            
        end %Trial
        
        C1=C1/Nt1;
        C2=C2/Nt2;
        
        %=== Try Alpha:
        for al=1:length(alphaSet) %try alphas
            clear W;clear W2; clear W1; clear L;
            alpha=alphaSet(al);
            
            [W1,L]=eig(C1,(C2+alpha*eye(size(C2))));
            [~,indx]=sort(diag(L),'descend');
            W1=W1(:,indx);
            clear L;
            
            [W2,L]=eig(C2,(C1+alpha*eye(size(C1))));
            [~,indx]=sort(diag(L),'descend');
            W2=W2(:,indx);
            W=[W1(:,1:m),W2(:,1:m)];
            
            %=== Training 1
            Ftrain1=nan(size(Train1,1),size(W,2));
            clear x; clear y; clear f;
            for Trial = 1:size(Train1,1)
                x=squeeze(Train1(Trial,:,:));
                x=x*W;
                f=var(x);
                Ftrain1(Trial,:)=f;
            end %Trial
            
            %=== Training 1
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
            PERF(K,al)= mean(pred==L_Test);
            
        end %al
        
    end %K
    
    performance=nanmean(PERF);
    [Perf_Sort,I_Sort]=sort(performance);
    Max_Perf = Perf_Sort(end);
    I = I_Sort(end);
    alpha=alphaSet(I);
else %find reqularized parameter alpha
    alpha=alphaSet;
end


%=== final W:
Nt1 = size(Class1,1);
Nt2 = size(Class2,1);

C1=0; C2=0;
for Trial = 1:max(Nt1,Nt2)
    
    if Trial<=Nt1
        x=squeeze(Class1(Trial,:,:));
        x = x- mean(x,1);
        c=x'*x/trace(x'*x);
        C1=C1+c;
    end %if
    
    if Trial<=Nt2
        x=squeeze(Class2(Trial,:,:));
        x = x- mean(x,1);
        c=x'*x/trace(x'*x);
        C2=C2+c;
    end %if
    
end %Trial

C1=C1/Nt1;
C2=C2/Nt2;

clear W1; clear W2; clear W;
[W1,L]=eig(C1,(C2+alpha*eye(size(C2))));
[~,indx]=sort(diag(L),'descend');
W1=W1(:,indx);
[W2,L]=eig(C2,(C1+alpha*eye(size(C1))));
[~,indx]=sort(diag(L),'descend');
W2=W2(:,indx);

WFinal=[W1(:,1:m),W2(:,1:m)];


end %function