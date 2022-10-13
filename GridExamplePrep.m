function [NodeInfo,NodeRegion,EdgeInfo,FlowClusterPre,FlowInfo,MinPath,PriorPairs,ExpCorrList]=GridExamplePrep(nNode)
%% NodeInfo [NodeID,xcoord,ycoord,P,A (for 3 scenarios)]
len=sqrt(nNode);
NodeInfo=(1:nNode)';
m=0; n=len;
for i=1:nNode
    m=m+1;
    NodeInfo(i,2:3)=[m,n];
    if m==len
        m=0; n=n-1;
    end
end

%% NodeRegion

NodeRegion=(1:nNode)';
NodeRegion(:,2)=[1;1;1;2;2;1;1;3;2;2;1;3;3;3;2;4;4;3;2;2;4;4;4;4;4];

%% EdgeInfo
nEdge=(len-1)*len*2;
EdgeInfo=[(1:nEdge);zeros(3,nEdge)];
k=1;
for m=1:nNode-1
    if rem(m,5)>0
        if (m+len)<=nNode
            EdgeInfo(2:3,k:k+1)=[m,m;m+1,m+len];
            k=k+2;
        else
            EdgeInfo(2:3,k)=[m;m+1];
            k=k+1;
        end
    else
        if (m+len)<=nNode
            EdgeInfo(2,k)=m;
            EdgeInfo(3,k)=m+len;
            k=k+1;
        end
    end
end
for i=1:nEdge
    EdgeInfo(4,i)=round(rand*10+5,2);
end

%% FlowClusterPre
nNode=25;
NodeRegion=(1:nNode)';
NodeRegion(:,2)=[1;1;1;2;2;1;1;3;2;2;1;3;3;3;2;4;4;3;2;2;4;4;4;4;4];
nFlowPre=nNode*(nNode-1)/2;
m=0;
FlowClusterPre=(1:nFlowPre)';
for i=1:nNode-1
    for j=i+1:nNode
        m=m+1;
        FlowClusterPre(m,2:3)=[i,j];
        FlowClusterPre(m,4:5)=[NodeRegion(NodeRegion(:,1)==i,2),NodeRegion(NodeRegion(:,1)==j,2)];
    end
end

for i=1:nFlowPre
    if FlowClusterPre(i,4)==FlowClusterPre(i,5)
        FlowClusterPre(i,6)=FlowClusterPre(i,4);
    else
        FlowClusterPre(i,6)=FlowClusterPre(i,4)+FlowClusterPre(i,5)+2;
        if FlowClusterPre(i,4)>1 && FlowClusterPre(i,5)>1
            FlowClusterPre(i,6)=FlowClusterPre(i,6)+1;
        end
    end
end
% X=FlowClusterPre(:,4:6); X=sortrows(X,'ascend');

%% FlowInfo (True mean, true std.dev)
FlowInfo=zeros(nFlowPre,9);
%% MinPath
Nodei=EdgeInfo(2,:);
Nodej=EdgeInfo(3,:);
Tij=EdgeInfo(4,:);

G=graph(Nodei,Nodej,Tij);	% graph with directed edges
plot(G, 'EdgeLabel', G.Edges.Weight);    % plot network with LinkID and c_ij

MinPath=struct; % Prepare struct for minimum paths
for i=1:nNode
    [MinPath(i).Paths,MinPath(i).T_ij] = shortestpathtree(G,i,'all','OutputForm','cell');
    %.Paths: minimum paths from i to other nodes
    %.T_ij: path cost from i to other nodes
    % LINK-OD INCIDENCE MATRIX & LINK ON MINIMUM PATHS
    for j=1:nNode
%         ODID=FlowInfo(FlowInfo(:,2)==i&FlowInfo(:,3)==j,1);
        len=size(MinPath(i).Paths{j,1},2);    % # of nodes to connect i and j
        if len>1 % origin and destination are different
            % IDENTIFY LINKS USED BY EACH MINIMUM PATH BRANCH
            MinPath(i).Links{j,1}=0;   % prepare array for included link IDs
            for k=1:len-1
                m=EdgeInfo(1,EdgeInfo(2,:)==MinPath(i).Paths{j,1}(k)&EdgeInfo(3,:)==MinPath(i).Paths{j,1}(k+1));
                if size(m,2)==0
                    m=EdgeInfo(1,EdgeInfo(2,:)==MinPath(i).Paths{j,1}(k+1)&EdgeInfo(3,:)==MinPath(i).Paths{j,1}(k));
                end    
                MinPath(i).Links{j,1}(k)=m;
            end

            % IDENTIFY OD PAIRS COVERED BY EACH MINIMUM PATH BRANCH
            MinPath(i).ODs{j,1}=0;   % prepare array for covered OD IDs
            col_ind=1;
            for k=1:len-1
                for m=k+1:len
                    if MinPath(i).Paths{j,1}(k)<MinPath(i).Paths{j,1}(m)
                        ODID2=FlowInfo(FlowInfo(:,2)==MinPath(i).Paths{j,1}(k)&FlowInfo(:,3)==MinPath(i).Paths{j,1}(m),1);
                    else
                        ODID2=FlowInfo(FlowInfo(:,2)==MinPath(i).Paths{j,1}(m)&FlowInfo(:,3)==MinPath(i).Paths{j,1}(k),1);
                    end
                    MinPath(i).ODs{j,1}(col_ind)=ODID2;
                    col_ind=col_ind+1;
                end
            end
        end
    end
end

%% True covariance matrix
% Rearrange FlowClusterPre
CorrPairCluster=struct;
ClusterN=unique(FlowClusterPre(:,6))';
ClusterN(2,:)=0;
m=0;
for i=1:nFlowPre
    if rand<0.7
        cluster=FlowClusterPre(i,6);
        ClusterN(2,cluster)=ClusterN(2,cluster)+1;
        CorrPairCluster(cluster).A(ClusterN(2,cluster),:)=FlowInfo(FlowInfo(:,1)==FlowClusterPre(i,1),:);
        m=m+1;
    end
end

TrCorrPairMean=[];  % prior information for correlated pairs
FlowCluster=[];     % cluster of pairs [ID, PairID, cluster]
for i=1:size(CorrPairCluster,2)
    TrCorrPairMean=[TrCorrPairMean;CorrPairCluster(i).A];
    B=CorrPairCluster(i).A(:,1);
    B(:,2)=i;
    FlowCluster=[FlowCluster;B];
end
nCorrFlow=size(FlowCluster,1);
FlowCluster=[(1:nCorrFlow)',FlowCluster];

UncorrFlow=setdiff(FlowClusterPre(:,1),FlowCluster(:,2));
nUncorr=length(UncorrFlow);
TrUncorr=[];
for i=1:nUncorr
    TrUncorr=[TrUncorr;FlowInfo(FlowInfo(:,1)==UncorrFlow(i,1),:)];
end

nBlocks=size(ClusterN,2); % # of available cluster combinations 
for i=1:nBlocks
    ClusterN(2,i)=sum(FlowCluster(:,3)==ClusterN(1,i));
end
BlockSize=ClusterN(2,:)';
RhoBlck=[0.4,0.5,0.3,0.3,0.2,0.3,0.5,0.2,0.3,0.4]; % base corr for clusters
delta=0.05;
eidim=4;
[TrueCorrMat]=CorrMatGen(nBlocks,BlockSize,RhoBlck,delta,eidim);
% nBlocks: # of blocks
% BlockSize: size of blocks (1-by-k array)
% RhoBlck: base constant correlation coefficients of blocks (1-by-k array)
% delta: base constant correlation coefficients of off-diagonal
% eidim: dimension of error space (>1)

% generate covariance matrix for scenarios
Truth=struct;
for scen=1:3
    Truth(scen).CorrMeanStd=TrCorrPairMean(:,[1,2*scen+2,2*scen+3]);
    Truth(scen).UncorrMeanStd=TrUncorr(:,[1,2*scen+2,2*scen+3]);
    D=zeros(nCorrFlow);
    for i=1:nCorrFlow
        D(i,i)=Truth(scen).CorrMeanStd(i,3); % assign std.dev to diagonal of D    
    end
    Truth(scen).CovMat=D*TrueCorrMat*D;   % derive covariance matrix
    PP=eig(Truth(scen).CovMat);
    sum(PP<0) 
end

%% PriorPairs
PriorPairs=zeros(nFlowPre,3);

%% ExpCorrList
% if equivalent
ExpCorrList=FlowCluster(:,2);
% if different
% ExpCorrList=SOMETHINGELSE
end