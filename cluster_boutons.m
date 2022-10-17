% Takes in a timeseries of noisy bouton data and clusters into putative
% axons. Input [DF] is a ROIs * frames matrix.

function [Clustered,IDs,Features,Models] = cluster_boutons(DF,varargin)
%% READY
Bin = 10;
NumFeatures = 40;
Ignore = false(size(DF,1),1);
Weighting = 'Skewness';
NumClustersPrecision = 10;

%% SET
for I=1:2:numel(varargin)
    eval([varargin{I} '= varargin{I+1};']);
end

%% GO
% get temporal features
[Features] = get_temporal_features(DF(:,~Ignore),'NumFeatures',NumFeatures,'BinSize',Bin);
    
% fit gmm
Options = statset('MaxIter',500,'TolFun',1e-6,'Display','final');

K = 1;
for NC = round(size(Features,1)/20):NumClustersPrecision:round(size(Features,1)/3) % range to test is reduced for speed
    TempModels{K} = fitgmdist(Features,NC,'RegularizationValue',1e-15,'CovarianceType','diagonal','Replicates',3,'Options',Options,'SharedCovariance',false,'ProbabilityTolerance',0,'Start','randSample');
    Criterion(K) = TempModels{K}.BIC;
    K = K +1;
end

[~,NumClustersID] = min(Criterion);
Model = TempModels{NumClustersID};

% cluster (1st round)
[IDs,~,Posteriors,~,~] = cluster(Model,Features);

[Clustered,Stray] = merge_clusters(DF,IDs,Posteriors,Weighting,1);

if ~any(Stray)
    StrayModel = [];
else
    % 2nd round
    StrayClusters = find(all((Clustered==0)'));
    
    StrayModel = fitgmdist(Features(Stray,:),length(StrayClusters),'RegularizationValue',1e-15,'CovarianceType','diagonal','Replicates',3,'Options',Options,'SharedCovariance',false,'ProbabilityTolerance',0,'Start','randSample');
    
    [StrayIDs,~,Posteriors,~,~] = cluster(StrayModel,Features(Stray,:));
    
    [StrayClustered] = merge_clusters(DF(Stray,:),StrayIDs,Posteriors,Weighting,2);
    
    % merge rounds
    Clustered(StrayClusters,:) = StrayClustered;
    StrayIDsID = unique(IDs(Stray));
    K = 1;
    NewStrayIDs = zeros(max(StrayIDs),1);
    for NewID = StrayIDsID'
        NewStrayIDs(StrayIDs==K) = NewID;
        K = K + 1;
    end
    IDs(Stray) = NewStrayIDs;
end

Models = {Model;StrayModel};
