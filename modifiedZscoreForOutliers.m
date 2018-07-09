% This program is to estimate the outliers and get rid of them.
% most of the outliers are located in the troughs so, each departure is set
% as zero. It is based on the Iglewicz and Hoaglin modified Z-score
% Comment 22/11/2011: the criterion changed to 1.2SD based on the
% comparisson between actual and retrieved noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input variables:
% [1] signalData = the ordinates of the BFP data
% [2] threshold = threshold to preserve the variance of the actual and split riples 

% Output variables: 
% [1] Position of the outliers
% [2] Filtered data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [outliersPosition,filteredSignalData] = modifiedZscoreForOutliers(signalData,threshold)
 % Initialy 3*theSD, then 1.20*theSD, then 1.10*SD
filteredSignalData = signalData;
theMean = mean(signalData);
theSD = std(signalData);
outliers = abs(signalData-theMean)>threshold*theSD; 

numberOfOutliers = sum(outliers);

outliersIndex = find(outliers>0);
outliersY = signalData(outliersIndex);
relatNegatOutliersIndex = find(outliersY<0);
absNegatOutliersIndex = outliersIndex(relatNegatOutliersIndex);
filteredSignalData(absNegatOutliersIndex) = -threshold*theSD;
outliersPosition = outliersIndex;




