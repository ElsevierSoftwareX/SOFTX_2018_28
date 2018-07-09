function [daughterWavelet, fourierFactor, coneOfInfl] = ...
    continuousWltBasis (motherWavelet, waveletScale, fourierFreq, wltParameter, sampleFreq)
% This function estimates the daughter of the mother wlt based on the
% parameters stated in the paper entitled A Practical Guide to Wavelet Analysis
% Compo & Torrence (1998)

% INPUTS
%     1. motherWavelet: the chosen mother wavelet,
%     2. waveletScale: the wavelet scale,
%     3. fourierFreq: a vector containing the Fourier Frequencies,  
%     4. wltParameter: the order of the derivative of the Gaussian (DOG) or
%        the central frequency of the Morlet wavelet
% 
% OUTPUTS
%     1. normalizationFactor: the factor estimated by appliying equation (6),
%     2. waveletScale: the scale of the wavelet transform estimated by applying eq. (9) and (10),
%     3. daughterWavelet: the conjugate of the Fourier transform of the mother Wavelet
    
%% Read General Data
stringMotherWlt = motherWavelet;
switch(stringMotherWlt)
    case {'DOG','MEXICANHAT'}, [daughterDOG, fourierFDOG, coneDOG] = ...
            DOGCalculator (waveletScale, fourierFreq, wltParameter, sampleFreq);
        daughterWavelet = daughterDOG;
        fourierFactor = fourierFDOG;
        coneOfInfl = coneDOG;
        
    case 'MORLET', [daughterMorlet, fourierFMorlet, coneMorlet] = ....
            morletCalculator (waveletScale, fourierFreq, wltParameter, sampleFreq);
        daughterWavelet = daughterMorlet;
        fourierFactor = fourierFMorlet;
        coneOfInfl = coneMorlet;
        
    otherwise disp('Not available Mother Wavelet');
end;
 
%% Calculator 
function [daughterDOG, fourierFDOG, coneDOG] = ...
    DOGCalculator (waveletScale, fourierFreq, wltParameter, sampleFreq)
normalizFactor = sqrt(2*pi/sampleFreq)/sqrt(gamma(wltParameter+0.5))*waveletScale^0.5;
expArgument = -(waveletScale.*fourierFreq).^2./2;
daughterDOG = -normalizFactor*(i^wltParameter)*...
    ((waveletScale.*fourierFreq).^wltParameter).*exp(expArgument);
fourierFDOG = 2*pi/sqrt(wltParameter+0.5);
coneDOG = fourierFDOG/sqrt(2); 
    
function [daughterMorlet, fourierFMorlet, coneMorlet] = ...
    morletCalculator (waveletScale, fourierFreq, wltParameter, sampleFreq)
centralFreq = wltParameter;
nFourierFreq=length(fourierFreq);
normalizFactor = sqrt(waveletScale*fourierFreq(2))*(pi^(-0.25))*sqrt(nFourierFreq);
expArgument = -(waveletScale.*fourierFreq-centralFreq).^2/2.*(fourierFreq>0);
daughterWavelet = normalizFactor*exp(expArgument);
daughterMorlet = daughterWavelet.*(fourierFreq>0);
fourierFMorlet = (4*pi)/(centralFreq+sqrt(2+centralFreq^2));
coneMorlet = fourierFMorlet/sqrt(2); 

