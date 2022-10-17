function [Features] = get_temporal_features(DFF,varargin)
%% READY
BinSize = 11; % ~half a second 
Method = 'ICA';
NumFeatures = 20;

%% SET
for I=1:2:numel(varargin)
    eval([varargin{I} '= varargin{I+1};']);
end

%% GO

% bin
NewTime = 1;
for Time = 1:BinSize:(size(DFF,2)-BinSize)
    Binned(:,NewTime) = nanmean(DFF(:,Time:Time+BinSize-1),2);
    NewTime = NewTime + 1;
end

% extract
if strcmp(Method,'ICA')
    Model = rica(Binned,NumFeatures,'IterationLimit',10000,'VerbosityLevel',0,'Standardize',true,'Lambda',0.0001);
    Features = transform(Model,Binned);
end
