%{
Knowledge Gradient with Independent Beliefs (KG)

notation for the following:
K: # of alternatives -> should it change per iteration
M: # of time-steps -> are we testing multiple time steps per iteration?
K x M stands for a matrix with K rows and M columns

This function takes in
mu:     true values for the mean (K x 1)
mu_0:   prior for the mean (K x 1)
sigma:  true values for std. dev (K x 1)
sigma_0 prior for std.dev (K x 1)
beta_W: measurement precision (1/lambda(x)) (K x 1)
M:      how many measurements will be made (scalar)

And returns
mu_est:     Final estimates for the means (K x 1)

OC:         Opportunity cost at each iteration (1 x M)
choices:    Alternatives picked at each iteration (1 x M)
mu_estALL:  Estimates at each iteration (K x M)
%}

% function [mu_est, OC, choices, mu_estALL, ObsOD, CovOD, PriorODMean]=Example_kg(mu,mu_0,M,ObsOD,CovOD,PriorODMean,CoveredOD,TrueODMean,TrueODCov)
% function [mu_est, OC, choices, ObsOD, PriorODMean, PriorODStdev]=Example_mab_transfer_20220720(mu,mu_0,M,ObsOD,PriorODMean,PriorODStdev,CoveredOD,ODSysNet,TrueODMean,TrueODStdev,TrueCov,CorFlow)
function [mu_est, OC, choices, ObsOD, PrUncorr, PrCorrMean,PrecisU,PrecisC]=Example_mab_transfer_20220801(mu,mu_0,M,ObsOD,CoveredOD,ODSysNet,PrUncorr,TrUncorr,PrCorrMean,TrCorrMean,TrCorrCov,PrecisU,PrecisC)
K=length(mu_0); %number of available choices

mu_est=mu_0;        % estimated mu from prior

OC=[];
choices=[];
% mu_estALL=[];
% nOD=size(PriorODMean,1);
nUncorr=size(PrUncorr,1);
nCorr=size(PrCorrMean,1);

% CREATE CovM, COVARIANCE MATRIX OF ACTIONS (CHOOSING A LINK), FROM CoveredOD
nCandLink=size(CoveredOD,2);

% % CREATE beta_W, ARRAY OF MEASUREMENT PRECISION (CHOICE LEVEL)
% beta_W=zeros(nCandLink,1);
% for i=1:nCandLink
%     beta_W(i,1)=1/sigma_0(i,1)^2;
% end
% beta=beta_W;

% RUN INITIAL nCandLink TRIALS TO TAKE AT LEAST ONE OBSERVATIONS FROM EACH
% ObsMAB=zeros(nCandLink,2);
ObsMAB=[mu_est,ones(nCandLink,1)];
% nCoverOD=zeros(nCandLink,1);
% 
% for k=1:nCandLink
%     nObs=size(ObsOD,2);
%     % observe result from takine each link
%     nCoverOD(k)=size(CoveredOD{k},2);
%     CovExtract=zeros(nCoverOD(k));
%     MeanExtract=zeros(nCoverOD(k),1);
%     for i=1:nCoverOD(k) % number of OD pairs covered by the best
%         subOD1=CoveredOD{k}(i);
%         for j=1:nCoverOD(k)
%             subOD2=CoveredOD{k}(j);
%             CovExtract(i,j)=TrueODCov(subOD1,subOD2);
%         end
%         MeanExtract(i,1)=TrueODMean(CoveredOD{k}(i),1);
%     end
%     W=mvnrnd(MeanExtract',CovExtract);
%     % update ObsOD
%     for i=1:nCoverOD(k)
%         ObsOD{CoveredOD{k}(i),nObs+1}=W(i);
%     end
%     ObsMAB(k,:)=[sum(W),1];
% end

% DETERMINE WHICH ACTION TO TAKE ACCORDING TO UCB ALGORITHM

for k=1:M % explore after nCandLink observation
    ArmValue=zeros(K,1);
    for iter1=1:K
        ArmValue(iter1)=ObsMAB(iter1,1)+sqrt(2*log(k)/ObsMAB(iter1,2));
    end
    
    [kBest,x]=max(ArmValue);
    
    %Plogy is the log values of KG for alternatives
    Plogy=[];
    Py=[];
%     KGAlts=zeros(K,8);
    
%     for iter1=1:K % fill 1st four columns
%         KGAlts(K,1)=mu_est(K,1);    % mean value
%         KGAlts(K,2)=beta(K,1);      % precision_n
%         KGAlts(K,3)=KGAlts(K,2)+beta_W(K,1);    % precision_n+1 
%         KGAlts(K,4)=sqrt(1/KGAlts(K,2)-1/KGAlts(K,3));  % change in variance
%     end
%     % 5th column (max value except for itself)
%     mu_est_sort=[(1:K)',mu_est];
%     [~,maxalt]=max(mu_est);
%     mu_est_sort=sortrows(mu_est_sort,2,'descend');
%     for iter1=1:K
%         if iter1==maxalt
%             KGAlts(K,5)=mu_est_sort(2,2);
%         else
%             KGAlts(K,5)=mu_est_sort(1,2);
%         end
%     end
%     
%     % 6-8th column
%     for iter1=1:K
%         KGAlts(K,6)=-abs((KGAlts(K,1)-KGAlts(K,5))/KGAlts(K,4)); % normalized influence of decision
%         KGAlts(K,7)=KGAlts(K,6)*normcdf(KGAlts(K,6))+normpdf(KGAlts(K,6)); % function of zeta
%         KGAlts(K,8)=KGAlts(K,4)*KGAlts(K,7);    % knowledge gradient
%     end
        
%         a=mu_est';  % choice-level mean
%         b=CovM(iter1,:)/sqrt(1/beta_W(iter1) + CovM(iter1,iter1));  % choice-level variation
%         if isreal(b)~=1 || size(a,2)~=size(b,2)
%             iter1
%         end
%         [KG,LogKG] = LogEmaxAffine(a,b);
%         
%         [Plogy]=[Plogy, LogKG];
%         [Py]=[Py,KG];
%     end
    
%     [~,x]=max(KGAlts(:,8));
    
    %max_value is the best estimated value of the KG 
    %x is the argument that produces max_value
    
    % ORIGINAL CODE
%     %observe the outcome of the decision
%     %W_k=mu_k+Z*SigmaW_k where SigmaW is standard deviation of the
%     %error for each observation
%     W_k=mu(x)+randn(1)*1./sqrt(beta_W(x)); % SHOULD BE DIFFERENT AS WE OBSERVE ALL AVAILABLE ODS

    %{
    HOW CAN WE UPDATE mu AND covM REGARDING DIFFERENT OBSERVATIONS OF ODS?
    I THINK WE HAVE TO DRAW RANDOM NUMBERS FROM EACH OD FLOW DISTRIBUTION
    AND UPDATE CovOD 1ST AND CovM 
    %}
    
    % RECEIVE OBSERVED OD FLOWS (RANDOM NUMBER REGARDING NORMAL DISTRIBUTION)
    % SINCE OUR PROBLEM INVOLVES DIFFERENT LEVEL OF FLOW, THE FUNCTION
    % SHOULD UPDATE THE INFORMATION ABOUT LOWER LEVEL (OD). UPDATE "ObsOD"
    % AND "PriorODMean".
    nObs=size(ObsOD,2);

    TotalOD = unique([CoveredOD{x},ODSysNet]);
    nCoverOD=size(TotalOD,2);
    
    MeanU=[];   p=0;
    StdevU=[];  
    InvolvedFlowU=[];
    InvolvedFlowC=[]; q=0;
    
    coveredODInd=zeros(nCoverOD,1);
    for i=1:nCoverOD
        if ismember(TotalOD(i),TrUncorr(:,1))
            p=p+1;
            InvolvedFlowU(p)=TotalOD(i);
            MeanU(p,1)=TrUncorr(TrUncorr(:,1)==TotalOD(i),2);
            StdevU(p,1)=TrUncorr(TrUncorr(:,1)==TotalOD(i),3);
        else
            q=q+1;
            InvolvedFlowC(q)=TotalOD(i);
        end
    end
    WU=normrnd(MeanU,StdevU);
    
    % Replace observation in W if correlated
    InvFl=size(InvolvedFlowC,2);
    MeanC=zeros(1,InvFl);
    CovC=zeros(InvFl);
    if InvFl>0
        for m=1:InvFl
%             MeanC(1,m)=TrueODMean(InvolvedFlow(m),1);
            MeanC(m)=TrCorrMean(TrCorrMean(:,1)==InvolvedFlowC(m),2);
            for n=m:InvFl
                % BRING COVARIANCES AND FORM A TEMPORARY COVARIANCE MATRIX
                CovC(m,n)=TrCorrCov(TrCorrMean(:,1)==InvolvedFlowC(m),TrCorrMean(:,1)==InvolvedFlowC(n));
                CovC(n,m)=CovC(m,n);
            end
        end
        WC=mvnrnd(MeanC,CovC);
    end
%     for m=1:InvFl
%         W(CoveredOD{x}==InvolvedFlow(m))=W2(m);
%     end
    
    % update ObsOD
    for i=1:p
        ObsOD{InvolvedFlowU(i),nObs+1}=WU(i);        
    end
    for i=1:q
        ObsOD{InvolvedFlowC(i),nObs+1}=WC(i);
    end
%     for i=1:nCoverOD
% %         ObsOD{CoveredOD{x}(i),nObs+1}=W(i);
%         ObsOD{TotalOD(i),nObs+1}=W(i);
%     end
    w_k=sum(WU)+sum(WC);
%     W_k=sum(W); % Observed total flow of choosing x
%     W2=W*coveredODInd;
%     W_k=sum(W2);
%     mu_est(x)=(mu_est(x)*(beta(x))+W_k*(beta_W(x)))/(beta(x)+beta_W(x)); % update mean
%     beta(x)=beta(x)+beta_W(x); % update precision


    ObsMAB(x,1)=(ObsMAB(x,1)*ObsMAB(x,2)+w_k)/(ObsMAB(x,2)+1);
    ObsMAB(x,2)=ObsMAB(x,2)+1;
    
%     e_x=zeros(K,1);       % prepare unit vector
%     e_x(x)=1;
%     addscalar = (W_k - mu_est(x))/(1/beta_W(x) + CovM(x,x));
%     mu_est=mu_est + addscalar*CovM*e_x;
%     CovM = CovM - (CovM*(e_x*e_x')*CovM)/((1/beta_W(x)) + CovM(x,x));
%     for i=1:nCandLink
%         if CovM(i,i)<0
%             k,i
%         end
%     end
    
    % % ORIGINAL CODE
%     e_x=zeros(K,1);
%     e_x(x)=1;
%     
%     %updating equations for Normal-Normal model with covariance
%     addscalar = (W_k - mu_est(x))/(1/beta_W(x) + covM(x,x));
%     mu_est=mu_est + addscalar*covM*e_x;
%     covM = covM - (covM*e_x*e_x'*covM)/((1/beta_W(x)) + covM(x,x));
    
    %pick the best one to compare OC
    [~, max_choice]=max(ObsMAB(:,1));

    %calculate the opportunity cost
    o_cost=max(mu)-mu(max_choice);
    
    OC=[OC,o_cost]; %update the OC matrix
    choices=[choices, x]; %update the choice matrix
    
    mu_est=ObsMAB(:,1);
%     if nargout>3 %if more than three outputs were asked
%         mu_estALL=[mu_estALL,mu_est];
%     end
    
    
end

% UPDATE "PriorODMean" and "CovOD" FOR OD PAIRS FROM ACCUMULATED "ObsOD"

nInvUncorr=length(InvolvedFlowU);
nInvCorr=length(InvolvedFlowC);
nExpUncorr=size(PrUncorr,1);
nExpCorr=size(PrCorrMean,1);
for i=1:nInvUncorr
    t=find(PrUncorr(:,1)==InvolvedFlowU(i));
    if t>0 % true uncorrelated flow is in uncorrelated prior
        if PrUncorr(t,3)>0
            precisN=1/PrUncorr(t,3)^2;
            PrUncorr(t,2)=(precisN*PrUncorr(t,2)+PrecisU(t,1)*WU(i))/(precisN+PrecisU(t,1));
            PrUncorr(t,3)=sqrt(1/(precisN+PrecisU(t,1)));
            PrUncorr(t,4)=PrUncorr(t,4)+1;
        else
            if WU(i)>0
                PrUncorr(t,2)=WU(i);
                PrUncorr(t,3)=WU(i)*0.01;
                PrUncorr(t,4)=PrUncorr(t,4)+1;
                PrecisU(t,1)=1/(PrUncorr(t,3)^2);
            end
        end
    else % true uncorrelated flow is in correlated prior
        v=find(PrCorrMean(:,1)==InvolvedFlowU(i));
        if PrCorrMean(v,3)>0
            precisN=1/PrCorrMean(v,3)^2;
            PrCorrMean(v,2)=(precisN*PrCorrMean(v,2)+PrecisC(v,1)*WU(i))/(precisN+PrecisC(v,1));
            PrCorrMean(v,3)=sqrt(1/(precisN+PrecisC(v,1)));
            PrCorrMean(v,4)=PrCorrMean(v,4)+1;
        else
            if WU(i)>0
                PrCorrMean(v,2)=WU(i);
                PrCorrMean(v,3)=WU(i)*0.01;
                PrCorrMean(v,4)=PrCorrMean(v,4)+1;
                PrecisC(v,1)=1/(PrCorrMean(v,3)^2);
            end
        end
    end
end

for i=1:nInvCorr
    t=find(PrCorrMean(:,1)==InvolvedFlowC(i));
    if t>0 % true correlated flow is in correlated prior
        if PrCorrMean(t,3)>0
            precisN=1/PrCorrMean(t,3)^2;
            PrCorrMean(t,2)=(precisN*PrCorrMean(t,2)+PrecisC(t,1)*WC(i))/(precisN+PrecisC(t,1));
            PrCorrMean(t,3)=sqrt(1/(precisN+PrecisC(t,1)));
            PrCorrMean(t,4)=PrCorrMean(t,4)+1;
        else
            if WC(i)>0
                PrCorrMean(t,2)=WC(i);
                PrCorrMean(t,3)=WC(i)*0.01;
                PrCorrMean(t,4)=PrCorrMean(t,4)+1;
                PrecisC(t,1)=1/(PrCorrMean(t,3)^2);
            end
        end
    else % true correlated flow is not in correlated prior
        v=find(PrUncorr(:,1)==InvolvedFlowC(i));
        if PrUncorr(v,3)>0
            precisN=1/PrUncorr(v,3)^2;
            PrUncorr(v,2)=(precisN*PrUncorr(v,2)+PrecisU(v,1)*WC(i))/(precisN+PrecisU(v,1));
            PrUncorr(v,3)=sqrt(1/(precisN+PrecisU(v,1)));
            PrUncorr(v,4)=PrUncorr(v,4)+1;
        else
            if WC(i)>0
                PrUncorr(v,2)=WC(i);
                PrUncorr(v,3)=WC(i)*0.01;
                PrUncorr(v,4)=PrUncorr(v,4)+1;
                PrecisU(v,1)=1/(PrUncorr(v,3)^2);
            end
        end
    end
end
    
% PriorODMean=zeros(nOD,2); % [mean,# of observation]
% InterOD=zeros(nOD,(size(ObsOD,2)-3));
% for i=1:nOD
%     for j=1:(size(ObsOD,2)-3)
%         if size(ObsOD{i,3+j},1)>0
%             InterOD(i,j)=1;
%         end
%     end
% end

% Flows=cell(nOD,1);
% for i=1:nOD
%     for j=4:size(ObsOD,2)
%         Flows{i,1}=[Flows{i,1},ObsOD{i,j}];
%     end
%     PriorODMean(i,2)=size(Flows{i,1},2);
%     if PriorODMean(i,2)>0
%         PriorODMean(i,1)=mean(Flows{i,1});
%     end
%     if PriorODMean(i,2)>1
%         PriorODStdev(i,1)=std(Flows{i,1});
%     end
% end

% UPDATE "PriorODStdev" FOR OD PAIRS FROM ACCUMULATED "ObsOD"
% CovOD=ones(nOD)*10+eye(nOD)*990;       % PREPARE COV MATRIX
% for i=1:nOD
%     for j=1:nOD
%         X=[];
%         Y=[];
%         CommonCase=InterOD(i,:).*InterOD(j,:);  % check if OD i and OD j are observed simultaneously
%         if sum(CommonCase)>1
%             for k=1:(size(ObsOD,2)-3)
%                 if CommonCase(k)==1
%                     X=[X,ObsOD{i,3+k}];
%                     Y=[Y,ObsOD{j,3+k}];
%                 end
%             end
%             COV=cov(X,Y);
%             CovOD(i,j)=COV(2,1);
%             CovOD(j,i)=CovOD(i,j);
%         end
%     end
% end


end

% logy = LogEmaxAffine(a,b)
% Calculates log(Exp[max_x a_x + b_x Z]-max_x a_x), where Z is a standard
% normal random variable and a,b are 1xM input vectors.
function [y,logy, a,b,c] = LogEmaxAffine(a,b)
    if (any(isnan(a)) || any(isnan(b)))
        warning('a or b is NaN');
    end
    assert(all(isreal(a)));
    assert(all(isreal(b)));

    a = a';
    b = b';

    % Check that a and b are column vectors of the right size
    if (any(size(a) ~= size(b)))
        error('LogEmaxAffine: a and b must be column vectors of the same size');
    end
    
    [a,b] = AffineBreakpointsPrep(a,b);
    
    [c, keep] = AffineBreakpoints(a,b); 
    a = a(keep);
    b = b(keep);
    c = c([1,keep+1]);
    M = length(keep);
    assert(all(isreal(c)));

    % I need logbdiff=log(diff(b)).  I thought that the following code would be
    % more numerically stable, able for example to distinguish cases like 
    % logb = [-25 -.3] vs. logb = [-35 -.3], but it doesn't seem to be able to.
    % Indeed, in the debugging output that I have below, the difference was 0.
    %{
    logb = log(abs(b)); % If b is 0, this is -Inf.
    sgnb = sign(b); % If b is 0, this is 0.
    logbdiff = zeros(size(c(2:M)));
    for i=1:length(b)-1
	[logbdiff(i),logbdiffsgn] = LogPlusExpSigned(logb(i),sgnb(i),logb(i+1),-sgnb(i+1));
	%assert(logbdiffsgn>=0);  % The b are distinct, so bdiff(i) can't be 0.
    end
    disp(sprintf('log(b)=%s log(diff(b))=%g logbdiff=%g difference=%g',mat2str(log(b)),log(diff(b)),logbdiff,log(diff(b))-logbdiff));
    %}
    logbdiff = log(diff(b))';  
    
    if M==1
        logy=log(a);
    elseif M>=2
        logy = LogSumExp(logbdiff+LogEI(-abs(c(2:M))));
    end

    logy=real(logy);
    y=exp(logy);
end
% Prepares vectors for passing to AffineEmaxBreakpoints, changing their
% order and removing elements with duplicate slope.

function [a,b] = AffineBreakpointsPrep(a,b)
    % Make sure a and b are column vectors.
    rows = size(a); if (rows == 1), a=a'; end
    rows = size(b); if (rows == 1), b=b'; end
    
    % 11/29/2008 PF: Experimental preprocessing step, which I hope will remove
    % a large number of the entries.
    [b1, i1] = min(b); % [a1,b1] is best at z=-infinity
    [a2, i2] = max(a); % [a2,b2] is best at z=0
    [b3, i3] = max(b); % [a3,b3] is best at z=+infinity
    a1 = a(i1);
    b2 = b(i2);
    a3 = a(i3);
    cleft = (a - a1)./(b1 - b); % intersection with leftmost line. 
    cright = (a - a3)./(b3 - b); % intersection with rightmost line.
    c2left = (a2 - a1)./(b1 - b2); % intersection with leftmost line. 
    c2right = (a2 - a3)./(b3 - b2); % intersection with rightmost line.
    keep = find(b==b1 | b==b3 | cleft <= c2left | cright >= c2right);
    %disp(sprintf('Preprocessing cut %d of %d entries', length(a)-length(keep), length(a)));
    a = a(keep);
    b = b(keep);
    clear keep cleft cright
   
        
    % Form a matrix for which ba(x,1) is the slope b(x) and ba(x,2) is the
    % y-intercept a(x).  Sort this matrix in ascending order of slope, 
    % breaking ties in slope with the y-intercept.  
    ba = [b, a];
    ba = sortrows(ba,[1,2]);
    a = ba(:,2);
    b = ba(:,1);
    
    % Then, from each pair of indices with the b component equal, remove
    % the one with smaller a component.  This code works because the sort
    % above enforced the condition: if b(i) == b(i+1), then a(i) <= a(i+1).
    keep = [find(diff(b)); length(b)];
    % This previous line is equivalent to:
    % keep = [];
    % for i=[1:length(b)-1]
    %    if b(i)~=b(i+1)
    %        keep = [keep, i];
    %    end
    %end 
    %keep = [keep, length(b)];  % We always keep the last one.
    
    % Note that the elements of keep are in ascending order.
    % This makes it so that b(keep) is still sorted in ascending order.
    a = a(keep);
    b = b(keep);
end

% Inputs are two M-vectors, a and b.
% Requires that the b vector is sorted in increasing order.
% Also requires that the elements of b all be unique.
% This function is used in AffineEmax, and the preparation of generic
% vectors a and b to satisfy the input requirements of this function are
% shown there.
%
% The output is an (M+1)-vector c and a vector A ("A" is for accept).  Think of
% A as a set which is a subset of {1,...,M}.  This output has the property
% that, for any i in {1,...,M} and any real number z,
%   i \in argmax_j a_j + b_j z
% iff
%   i \in A and z \in [c(j+1),c(i+1)],
%   where j = sup {0,1,...,i-1} \cap A.
%
% A note about indexing:
% Since Matlab does not allow indexing from 0, but instead requires
% indexing from 1, what is called c_i in the paper is written in matlab as
% c(1+i).  This is because in the paper we reference c_0.  For the vectors
% a and b, however, we don't need to reference a_0 or b_0, so we reference
% a_i and b_i by a(i) and b(i) respectively, rather than a(i+1) or b(i+1).
% 
function [c,A] = AffineBreakpoints(a,b)
    % Preallocate for speed.  Instead of resizing the array A whenever we add
    % to it or delete from it, we keep it the maximal size, and keep a length
    % indicator Alen telling us how many of its entries are good.  When the
    % function ends, we remove the unused elements from A before passing
    % it.
    M = length(a);
    c = zeros(1,M+1);
    A = zeros(1,M);
    
    % Step 0
    i=0;
    c(1+i) = -inf;
    c(1+i+1) = +inf;
    A(1) = 1;
    Alen = 1;
    
    for i=[1:M-1]
        c(1+i+1) = +inf;
        while(1)
            j = A(Alen); % jindex = Alen
            c(1+j) = (a(j) - a(i+1))/(b(i+1)-b(j));
	    % The if statement below replaces these lines from version 2 of the
	    % function.
	    %    kindex = jindex-1 = Alen-1
            %    if kindex > 0 && c(1+j)<=c(1+A(kindex))
            if Alen > 1 && c(1+j)<c(1+A(Alen-1))
		Alen = Alen-1; % Remove last element j
                % continue in while(1) loop
            else
                break % quit while(1) loop
            end
        end
	A(Alen+1) = i+1;
	Alen = Alen + 1;
    end
    A = A(1:Alen);
end

% Returns the log of Exp[(s+Z)^+], where s is a constant and Z is a standard
% normal random variable.  For large negative arguments Exp[(s+Z)^+] function
% is close to 0.  For large positive arguments, the function is close to the
% argument.  For s large enough, s>-10, we use the formula
% Exp[(s+Z)^+] = s*normcdf(s) + normpdf(s).  For smaller s we use an asymptotic
% approximation based on Mill's ratio.  EI stands for "expected improvement",
% since Exp[(s+Z)^+] would be the log of the expected improvement by measuring
% an alternative with excess predictive mean s over the best other measured
% alternative, and predictive variance 0.
function logy = LogEI(s)

% Use the asymptotic approximation for these large negative s.  The
% approximation is derived via:
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-|s|normcdf(-|s|)/normpdf(s)]
% and noting that normcdf(-|s|)/normpdf(s) is the Mill's ratio at |s|, which is
% asymptotically approximated by |s|/(s^2+1) [Gordon 1941, also documented in
% Frazier,Powell,Dayanik 2009 on page 14].  This gives,
%   s*normcdf(s) + normpdf(s) = normpdf(s)*[1-s^2/(s^2+1)] = normpdf(s)/(s^2+1).

i=find(s<-10);
if (length(i)>0)
    logy(i) = LogNormPDF(s(i)) - log(s(i).^2 + 1);
end

% Use straightforward routines for s in the more numerically stable region.
i=find(s>=-10);
if (length(i)>0)
    logy(i) = log(s(i).*normcdf(s(i))+normpdf(s(i)));
end

assert(all(isreal(logy)));
end

% logy = LogNormPDF(z)
% Returns the log of the normal pdf evaluated at z.  z can be a vector or a scalar.
function logy = LogNormPDF(z)
	const = -.5*log(2*pi); % log of 1/sqrt(2pi).
	logy = const - z.^2/2;
end

% function y=LogSumExp(x)
% Computes log(sum(exp(x))) for a vector x, but in a numerically careful way.
function y=LogSumExp(x)
xmax = max(x);
y = xmax + log(sum(exp(x-xmax)));
end
