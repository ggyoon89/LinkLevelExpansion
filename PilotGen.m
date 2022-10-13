function [PrUncorrPair,PrCorrPairMean,PrCorrFlowCov,ObsOD,PilotStartEnd,PrecisU,PrecisC]=PilotGen(nPilot,nTimes,minpilotL,nNode,nOD,MinPath,TrueUncorrPair,TrueCorrPairMean,TrueCorrPairCov,FlowCluster,PriorPairs,ExpCorrList)
% function [PrUncorrPair,PrCorrPairMean,PrCorrFlowCov,ObsOD,PilotStartEnd,PrecisU,PrecisC]=PilotGen(nPilot,nTimes,minpilotL,nNode,nOD,MinPath,TrueUncorrPair,TrueCorrPairMean,TrueCorrPairCov,FlowCluster,PrUncorrPair,PrCorrPairMean,PilotStartEnd)
% % INPUTS
% nPilot: # of pilots
% nTimes: # of operations per pilot
% ODinfo: information of OD pairs [ID,i,j,d_ij,X_ij,std.dev(X_ij)]
% TrueODCov: true covariance matrix of OD flows

%% SYNTHESIZE COVARIANCE MATRIX WITH COVARIANCE OF OD FLOWS
% LET'S OBSERVE X_ij* FROM PILOTS BY PICKING A MINIMUM PATH ROUTE TO OPERATE
ObsOD=num2cell((1:nOD)');  % ID of observed OD flow [OD ID]
ObsOD{1,1+nPilot}=[];      % empty columns for observed OD flows
TUCF=TrueUncorrPair;
TCFMean=TrueCorrPairMean;
TCFCov=TrueCorrPairCov;
TrCorFlows=FlowCluster(:,2);
PrUncorrPair=[];
PrCorrPairMean=[];
% split PriorPairs
for i=1:nOD
    if ismember(PriorPairs(i,1),ExpCorrList)
        PrCorrPairMean=[PrCorrPairMean;PriorPairs(i,:)];
    else
        PrUncorrPair=[PrUncorrPair;PriorPairs(i,:)];
    end
end
nExpCorrPair=size(PrCorrPairMean,1);
nExpUncorrPair=size(PrUncorrPair,1);

for i=1:nPilot
%     j=PilotStartEnd(i,1); % start node of pilot route, j
%     k=PilotStartEnd(i,2); % end node of pilot route, k
    j=0;    % start node of pilot route, j
    k=0;    % end node of pilot route, k
    while j==k % when an intrazonal trip is chosen
        j=randi(nNode); % randomize j
        k=randi(nNode); % randomize k
        if j~=k % when an interzonal trip is chosen
            npathOD=size(MinPath(j).ODs{k,1},2);    % # of ODs covered by minimum path j to k
            if npathOD>=2*nchoosek(minpilotL,2) % choose path that consists of more than minpilotL nodes
                % OBSERVE OD FLOWS BASED ON COVARIANCE MATRIX
                MeanFlow=zeros(1,npathOD);   % array for true mean OD flows
                CovFlow=zeros(npathOD);      % matrix for true covariance
                ODSet=MinPath(j).ODs{k,1};  % array for flow ID replacement
                for m=1:npathOD
                    if ismember(ODSet(m),TrCorFlows) % true mean OD flow of m-th OD
                        MeanFlow(m)=TCFMean(TrCorFlows==ODSet(m),2);
                    else
                        MeanFlow(m)=TUCF(TUCF(:,1)==ODSet(m),2);
                    end
%                     MeanFlow(m)=ODinfo(IDFlow(m),4);   % true mean OD flow of m-th OD
                    for n=m:npathOD
                        if ismember(ODSet(m),TrCorFlows)&&ismember(ODSet(n),TrCorFlows)
                            mID=FlowCluster(TrCorFlows==ODSet(m),1);
                            nID=FlowCluster(TrCorFlows==ODSet(n),1);
                            CovFlow(m,n)=TCFCov(mID,nID);  % true covariance of m- and n-th OD
                            CovFlow(n,m)=CovFlow(m,n);  % symmetry
                        elseif m==n
                            CovFlow(m,n)=TUCF(TUCF(:,1)==ODSet(m),3)^2;
                        end
                    end
                end
                Obs=round(mvnrnd(MeanFlow,CovFlow,nTimes)); % array of observation of i-th pilot
                if sum(Obs<0)>0 % check if any negative flow is observed
                    i
                end
                for m=1:npathOD % input observed flows to ObsOD
                    ObsOD{ODSet(m),1+i}=Obs(:,m)';
                end
            else
                j=k; % repeat the loop 
            end
            PilotStartEnd(i,1:2)=[j,k];
        end
    end
end

% UPDATE PRIOR FROM OBSERVATION
% 0. Precision
% Assume as 1% of prior mean
PrecisU=zeros(nExpUncorrPair,1);
PrecisC=zeros(nExpCorrPair,1);
for i=1:nExpUncorrPair
    if PrUncorrPair(i,2)>0
        PrecisU(i,1)=1/(PrUncorrPair(i,2)*0.01)^2;
    else
        PrecisU(i,1)=1;
    end
end
for i=1:nExpCorrPair
    if PrCorrPairMean(i,2)>0
        PrecisC(i,1)=1/(PrCorrPairMean(i,2)*0.01)^2;
    else
        PrecisC(i,1)=1;
    end
end

% I. Uncorrelated flow
% Update mean and build std.dev of uncorrelated flows
% Assume we don't know both true mean and variance are unknown
% Conjugate prior distribution: Normal-Gamma
% 1. Update mean (priors are based on single observation)
% 2. Build std.dev

for i=1:nExpUncorrPair
    K=cell2mat(ObsOD(PrUncorrPair(i,1),2:end));
    k=length(K);
    if k>0
        PrUncorrPair(i,2)=(PrUncorrPair(i,2)+sum(K))/(1+k);
        PrUncorrPair(i,3)=std([PrUncorrPair(i,2),K]);
        PrUncorrPair(i,4)=1+k;
        PrecisU(i,1)=1/(PrUncorrPair(i,2)*0.01)^2;
    end
end

% II. Correlated flow
% Update mean and build covariance matrix of correlated flows
% Since we don't have prior covariance, just simply update mean.

for i=1:nExpCorrPair
    K=cell2mat(ObsOD(PrCorrPairMean(i,1),2:end));
    k=length(K);
    if k>0
        PrCorrPairMean(i,2)=(PrCorrPairMean(i,2)+sum(K))/(1+k);
        % Column3 is std.dev used for true covariance
        PrCorrPairMean(i,4)=1+k;
        PrecisC(i,1)=1/(PrCorrPairMean(i,2)*0.01)^2;
    end
end

% Prior covariance is made from observed pilots
% Assume measurement error of 1% of mean: apply to unobserved 
% Incidence matrix among correlated flows in pilot
IncCorrFlow=zeros(nExpCorrPair,nPilot);  % ODs observed in pilots
for i=1:nExpCorrPair
    for j=1:nPilot
        if size(ObsOD{i,1+j},1)>0
            IncCorrFlow(i,j)=1; % i-th OD pair is served by j-th pilot
        end
    end
end

% Initiate COVARIANCE MATRIX OF OD FLOW (nCorrFlow-BY-nCorrFlow)
% Assume measurement error of 1% of mean
PrCorrFlowCov=zeros(nExpCorrPair);
% for i=1:nCorrPair
%     PrCorrFlowCov(i,i)=PrecisC(i,1);
% end
CovODInd=zeros(nExpCorrPair);    % prepare matrix indicating the existence of covariance
for i=1:nExpCorrPair
    for j=i:nExpCorrPair
        X=[];   % observation of OD flow i
        Y=[];   % observation of OD flow j
        CommonCase=IncCorrFlow(i,:).*IncCorrFlow(j,:);  % check if OD i and OD j are observed simultaneously
        if sum(CommonCase)>1
            for k=1:nPilot
                if CommonCase(k)==1 % both OD i and j are served at the same time
                    X=[X,ObsOD{i,1+k}];
                    Y=[Y,ObsOD{j,1+k}];
                end
            end
            COV=cov(X,Y);           % covariance matrix of X and Y (2-by-2 matrix)
            PrCorrFlowCov(i,j)=COV(2,1);    % covariance
            CovODInd(i,j)=1;        % indicator
            PrCorrFlowCov(j,i)=PrCorrFlowCov(i,j);  % symmetry
            CovODInd(j,i)=CovODInd(i,j);    % symmetry
        else
            PrCorrFlowCov(i,j)=0;
            PrCorrFlowCov(j,i)=0;
        end
    end
    if PrCorrFlowCov(i,i)==0
        PrCorrFlowCov(i,i)=1/sqrt(PrecisC(i,1));
    end
end

for i=1:nExpCorrPair
    PrCorrPairMean(i,3)=sqrt(PrCorrFlowCov(i,i)); % update prior std.dev from covariance matrix
end

% % Derive incidence matrix
% InterOD=zeros(nOD,nPilot);  % ODs observed in pilots
% for i=1:nOD
%     for j=1:nPilot
%         if size(ObsOD{i,1+j},1)>0
%             InterOD(i,j)=1; % i-th OD pair is served by j-th pilot
%         end
%     end
% end
% 
% % DERIVE COVARIANCE MATRIX OF OD FLOW (nOD-BY-nOD)
% % CovOD=ones(nOD)*10^(3)+eye(nOD)*9000;       % prepare initial belief
% % CovOD=eye(nOD)*100;
% CovOD=zeros(nOD);
% for i=1:nOD
%     CovOD(i,i)=(PriorODMean(i,1)*0.1)^2;
% end
% CovODInd=zeros(nOD);    % prepare matrix indicating the existence of covariance
% for i=1:nOD
%     for j=1:nOD
%         X=[];   % observation of OD flow i
%         Y=[];   % observation of OD flow j
%         CommonCase=InterOD(i,:).*InterOD(j,:);  % check if OD i and OD j are observed simultaneously
%         if sum(CommonCase)>1
%             for k=1:nPilot
%                 if CommonCase(k)==1 % both OD i and j are served at the same time
%                     X=[X,ObsOD{i,1+k}];
%                     Y=[Y,ObsOD{j,1+k}];
%                 end
%             end
%             COV=cov(X,Y);           % covariance matrix of X and Y (2-by-2 matrix)
%             CovOD(i,j)=COV(2,1);    % covariance
%             CovODInd(i,j)=1;        % indicator
%             CovOD(j,i)=CovOD(i,j);  % symmetry
%             CovODInd(j,i)=CovODInd(i,j);    % symmetry
%         end
%     end
% end
% 
% PriorODStdev=zeros(nOD,1);   % [mean,# of observation]
% for i=1:nOD
%     if CovOD(i,i)>0
%         PriorODStdev(i,1)=sqrt(CovOD(i,i)); % update prior mean from covariance matrix
%     end
% end