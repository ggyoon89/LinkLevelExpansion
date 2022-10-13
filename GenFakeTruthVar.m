function [TrueUncorrPair,TrueCorrPairMean,TrueCorrPairCovL,TrueCorrPairCovM,TrueCorrPairCovH,FlowCluster,PrUncorrPair,PrCorrPairMean]=GenFakeTruthVar(PUMAinfo,nPairSignif,NodeCluster)
%{
[OUTPUTS]
TrueUncorrFlow: fake true mean and std.dev of uncorrelated flows (order of FlowID)
TrueCorrFlowMean: fake true mean of corrrelated flows (order of paired flow level, both directions)
TrueCorrFlowCov: fake true covariance matrix of correlated flows
FlowCluster: clusters of flows (1~7)
PrUncorrFlow: mean and std.dev prior of uncorrelated flows
PrCorrFlowMean: mean prior of correlated flows (no covariance matrix at the beginning) + true std.dev

[GENERATE FAKE TRUTH FROM ACTUAL DATASET]
1. Accept actual RHTS result (flows) as our prior knowledge
2. Assume that "Prior knowledge may be originated from fake truth."
3. Randomly generate fake truths (mean, standard deviation, and covariance)
 3-1. Randomly generate "artificial" standard deviation vector for other
      flows. They can be considered as independent variables. 
 3-2. Randomly generate "artificial" mean vector for flows. For those from
      3-1, use "mvnrnd." For others, use "normrnd."
 3-3. Randomly generate "artificial" covariance matrix for a subset of 
      flows. This does not have to be equivalent to the "significant flows"
      designated by operator.
4. The result consists of 4 priors.
 4-1. Mean of correlated flows
 4-2. Covariance of correlated flows
 4-3. Mean of other flows
 4-4. Standard deviation of other flows
%}

P=PUMAinfo;
NC=NodeCluster;
Pairs=sortrows(P,4,'descend'); % sort (994 out of 1485 have flow)
% nPairwithFlow=sum(PairedOD(:,8)>0);
% nPairSignif=round(nPairwithFlow/2); % Assume top 50% are correlated
% nPairSignif=500;
CorrelatedPairs=Pairs(1:nPairSignif,1)';

%% Generate std.dev of flows
Mu0=P(:,4); % set actual RHTS result (flows) as our mean prior
nPair=length(Mu0);
for i=1:nPair
    if P(i,4)>0
        % Range of std.dev differs by the number of digits
        maxdev=((6-log10(P(i,4)))*2+3)/100; mindev=0.02; % low (2-(3-13%))
        P(i,5)=P(i,4)*(mindev+(maxdev-mindev)*rand);
        maxdev=((6-log10(P(i,4)))*2+13)/100; mindev=0.10; % mid (10-(13-23%))
        P(i,6)=P(i,4)*(mindev+(maxdev-mindev)*rand);
        maxdev=((6-log10(P(i,4)))*2+23)/100; mindev=0.20; % high (20-(23-33%))
        P(i,7)=P(i,4)*(mindev+(maxdev-mindev)*rand);
    end
end

% Separate flows 
PrCorrPairMean=zeros(1,5);  % [FlowID, mean, std.dev (l,m,h)]
PrUncorrPair=zeros(1,5);    % [FlowID, mean, std.dev (l,m,h)]
m=0;n=0;
for i=1:nPair
    if ismember(i,CorrelatedPairs)
        m=m+1;
        PrCorrPairMean(m,:)=P(i,[1,4:7]);
    else
        n=n+1;
        PrUncorrPair(n,:)=P(i,[1,4:7]);
    end
end
nCorrPair=m; nUncorrPair=n;
%% Generate fake truth (mean, std.dev, covariance)
% Uncorrelated flow: norminv
TrueUncorrPair=PrUncorrPair; % save true information
PrUncorrPair(:,3)=PrUncorrPair(:,2)*0.05; % assume prior std.dev as 5% of prior mean
TrueUncorrPair(:,2)=norminv(rand(nUncorrPair,1),PrUncorrPair(:,2),TrueUncorrPair(:,3)); % generate true mean
for i=1:nUncorrPair
    if TrueUncorrPair(i,3)==0
        TrueUncorrPair(i,2)=0;
    end
    if TrueUncorrPair(i,2)<0
        TrueUncorrPair(i,2)=norminv(rand,PrUncorrPair(i,2),TrueUncorrPair(:,3));
    end
end
% TrueUncorrFlow: [FlowID, mean, std.dev (l,m,h)] in FlowID order
Q=PrCorrPairMean;
%% Correlated flow
% rearrange PrCorrFlowMean
% before: order of paired OD flows (1a,1b),(2a,2b),...
% after: two sequence of ordered flows (1a,2a,...),(1b,2b,...)

% PrCorrFlowMean(:,:)=0;
NCluster=zeros(1,7);
CorrPairCluster=struct;
for i=1:nCorrPair
    ClusterPair=[NC(P(P(:,1)==Q(i,1),2),2),NC(P(P(:,1)==Q(i,1),3),2)];
    if ClusterPair(1)==ClusterPair(2)
        cluster=ClusterPair(1); % within cluster (1/2/3/4)
    else
        b=ClusterPair(1)*ClusterPair(2);
        if b==3 % 1 & 3
            cluster=5;
        elseif b==4 % 1 & 4
            cluster=6;
        elseif b==12 % 3 & 4
            cluster=7;
        else % 1 & 2 / 2 & 3 / 2 & 4
            cluster=2; % aggregate to Cluster 2
        end
    end
    NCluster(cluster)=NCluster(cluster)+1;
    CorrPairCluster(cluster).A(NCluster(cluster),:)=Q(Q(:,1)==CorrelatedPairs(i),:);
%     CorrPairCluster(cluster+7).A(NCluster(cluster),:)=Q(Q(:,1)==CorrelatedPairs(2*i),:);
%     PrCorrFlowMean(i,:)=Q(Q(:,1)==CorrelatedFlows(2*i-1),:);
%     PrCorrFlowMean(i+nCorrFlow/2,:)=Q(Q(:,1)==CorrelatedFlows(2*i),:);
end

%% reassemble PrCorrFlowMean
PrCorrPairMean=[];
FlowCluster=[];
for i=1:size(CorrPairCluster,2)
    PrCorrPairMean=[PrCorrPairMean;CorrPairCluster(i).A];
%     if i>7
%         j=i-7;
%     else
%         j=i;
%     end
    B=CorrPairCluster(i).A(:,1);
    B(:,2)=i;
    FlowCluster=[FlowCluster;B];
end
FlowCluster=[(1:size(FlowCluster,1))',FlowCluster];
%% generate correlation matrix
nBlocks=7; % # of availble cluster combinations 
BlockSize=NCluster;
% BlockSize=[100,100,100,100,97];
RhoBlck=[0.4,0.5,0.3,0.3,0.2,0.3,0.5];
% RhoBlck=[0.4,0.2,0.3,0.5,0.6];
delta=0.05;
eidim=6;
[TrueCorrMat]=CorrMatGen(nBlocks,BlockSize,RhoBlck,delta,eidim);
% nBlocks: # of blocks
% BlockSize: size of blocks (1-by-k array)
% RhoBlck: base constant correlation coefficients of blocks (1-by-k array)
% delta: base constant correlation coefficients of off-diagonal
% eidim: dimension of error space (>1)

% generate covariance matrix
DL=zeros(nCorrPair);
DM=zeros(nCorrPair);
DH=zeros(nCorrPair);
for i=1:nCorrPair
    DL(i,i)=PrCorrPairMean(i,3); % assign std.dev to diagonal of D  
    DM(i,i)=PrCorrPairMean(i,4);
    DH(i,i)=PrCorrPairMean(i,5);
end
TrueCorrPairCovL=DL*TrueCorrMat*DL;   % derive covariance matrix
TrueCorrPairCovM=DM*TrueCorrMat*DM;
TrueCorrPairCovH=DH*TrueCorrMat*DH;
PP=eig(TrueCorrPairCovL);
sum(PP<0) 
PP=eig(TrueCorrPairCovM);
sum(PP<0) 
PP=eig(TrueCorrPairCovH);
sum(PP<0) 

% generate truth 
TrueCorrPairMean=[PrCorrPairMean(:,1),mvnrnd(PrCorrPairMean(:,2),TrueCorrPairCovL)',PrCorrPairMean(:,3:5)];
PrCorrPairMean(:,3)=PrCorrPairMean(:,2)*0.05;
PrCorrPairMean(:,4)=PrCorrPairMean(:,2)*0.05;
PrCorrPairMean(:,5)=PrCorrPairMean(:,2)*0.05;

end