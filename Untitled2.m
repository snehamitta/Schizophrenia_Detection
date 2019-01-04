function signal_filt = Untitled2(signal,SamplePeriod,LowCutoff,HighCutoff,AddMeanBack)
% 
% function signal_filt = bandpass_ideal_rect_wdw(...
%     signal,SamplePeriod,LowCutoff,HighCutoff,AddMeanBack)
% 
% Ideal rectangular window for bandpass filtering. Useful for resting state
% studies.
% 
% 
%INPUT
%-----
% - SIGNAL      : 2D data matrix (nTimePoints x nTimeSeries)
% - SAMPLEPERIOD: sampling period [s]
% - LOWCUTOFF   : high pass, low cutoff of the band [Hz]
%   (use 0 for low-pass filter)
% - HIGHCUTOFF  : low pass, high cutoff of the band [Hz]
%   (use 0 for high-pass filter)
% - ADDMEANBACK : 1 (add the mean back after filtering) or 0 (do not add)
% 
% 
%OUTPUT
%------
% - SIGNAL_FILT: filtered signal
% 
% 
%EXAMPLE
%-------
% T  = .5; % [s] => sampling freq. = 1/T => Nyquist freq. = 1/(2T)
% t  = (0:T:2*pi*10)';
% f1 = (1/2/T)*.030;
% f2 = (1/2/T)*.050;
% f3 = (1/2/T)*.075;
% x1 = 10*sin(2*pi*f1*t);
% x2 = 4*sin(2*pi*f2*t);
% x3 = 25*cos(2*pi*f3*t);
% x = x1 + x2 + x3;
% xLP = bandpass_ideal_rect_wdw(x,T,0,(f2+f3)/2,1);
% figure, plot(t,x,t,xLP,t,x1+x2)
% title('Low-pass filter'), legend('Original', 'Filtered', 'f1 & f2')
% xHP = bandpass_ideal_rect_wdw(x,T,(f1+f2)/2,0,1);
% figure, plot(t,x,t,xHP,t,x2+x3)
% title('High-pass filter'), legend('Original', 'Filtered', 'f2 & f3')
% xBP = bandpass_ideal_rect_wdw(x,T,(f1+f2)/2,(f2+f3)/2,1);
% figure, plot(t,x,t,xBP,t,x2)
% title('Band-pass filter'), legend('Original', 'Filtered', 'f2')
% 
% 
% Based on REST_BANDPASS and REST_IDEALFILTER

% - CUTNUMBER: cut the data into pieces if small RAM memory, e.g., 4GB is
%   available on PC. It can be set to 1 on server with big memory (e.g.,
%   50GB). Default: 10

% Guilherme Coco Beltramini - 2012-Oct-23, 06:25 pm
% 
% Changelog
%--------------------------------------------------------------------------
% 2012-Oct-24, 10:37 am: "cutnumber" removed. Only useful when the number
%   of time series is huge (fMRI studies). Set to 10 for fMRI studies.

cutnumber = 1;


% Remove the mean
%==========================================================================
if size(signal,1)==1
    signal = signal.';
end
%nTimePoints = size(signal,1);
%theMean     = mean(signal,1);
%theMean     = repmat(theMean,nTimePoints,1);
%signal      = signal - theMean;
theMean = sum(signal, 1)./size(signal, 1);
signal  = bsxfun(@minus, signal, theMean);


SegmentLength = ceil(size(signal,2)/cutnumber);
for iCut=1:cutnumber
    
    % Find the segment
    %-----------------
    if iCut~=cutnumber
        Segment = (iCut-1)*SegmentLength+1 : iCut*SegmentLength;
    else
        Segment = (iCut-1)*SegmentLength+1 : size(signal, 2);
    end
    
    % Filter the data
    %----------------
    signal(:,Segment) = ideal_filter(signal(:,Segment),SamplePeriod,[LowCutoff,HighCutoff]);
    
end


% Add the mean back after filter
%==========================================================================
if AddMeanBack
    %signal = signal + theMean;
    signal = bsxfun(@plus, signal, theMean);
end

signal_filt = signal;

end


function [Data_Filtered] = ideal_filter(Data,SamplePeriod,Band)
% 
% [Data_Filtered] = ideal_filter(Data,SamplePeriod,Band)
% 
% 
%INPUT
%-----
% - Data        : 2D data matrix (nDimTimePoints x nTimeSeries)
% - SamplePeriod: Sample period, i.e., 1/sample frequency
% - Band        : frequency for filtering (1x2 array):
%   [LowCutoff_HighPass HighCutoff_LowPass]: band-pass filtering
%   [0 HighCutoff_LowPass]: low-pass filtering
%   [LowCutoff_HighPass 0]: high-pass filtering
% 
% 
%OUTPUT
%------
% - Data_Filtered: data after filtering
% 

% Written by YAN Chao-Gan 120504 based on REST.
% The Nathan Kline Institute for Psychiatric Research, 140 Old Orangeburg Road, Orangeburg, NY 10962, USA
% Child Mind Institute, 445 Park Avenue, New York, NY 10022, USA
% The Phyllis Green and Randolph Cowen Institute for Pediatric Neuroscience, New York University Child Study Center, New York, NY 10016, USA
% ycg.yan@gmail.com


sampleFreq 	 = 1/SamplePeriod;
sampleLength = size(Data,1);
paddedLength = nextpow2_one35(sampleLength); %2^nextpow2(sampleLength);
LowCutoff_HighPass = Band(1);
HighCutoff_LowPass = Band(2);


% Get the frequency index
%==========================================================================
if (LowCutoff_HighPass >= sampleFreq/2) % All high stop
    idxLowCutoff_HighPass = paddedLength/2 + 1;
else % high pass, such as freq > 0.01 Hz
    idxLowCutoff_HighPass = ceil(LowCutoff_HighPass * paddedLength * SamplePeriod + 1);
end

if (HighCutoff_LowPass>=sampleFreq/2) || (HighCutoff_LowPass==0) % All low pass
    idxHighCutoff_LowPass = paddedLength/2 + 1;
else % low pass, such as freq < 0.08 Hz
    idxHighCutoff_LowPass = fix(HighCutoff_LowPass * paddedLength * SamplePeriod + 1);
end

FrequencyMask = zeros(paddedLength,1);
FrequencyMask(idxLowCutoff_HighPass:idxHighCutoff_LowPass,1) = 1;
FrequencyMask(paddedLength-idxLowCutoff_HighPass+2:-1:paddedLength-idxHighCutoff_LowPass+2,1) = 1;

FrequencySetZero_Index = FrequencyMask==0;


% FFT -> filter -> IFFT
%==========================================================================

% Remove the mean before zero padding
%------------------------------------
%Data = Data - repmat(mean(Data,1), size(Data,1), 1);
Data = bsxfun(@minus, Data, sum(Data, 1)./size(Data, 1));

% Pad with zero
%--------------
Data = [Data ; zeros(paddedLength-sampleLength, size(Data, 2))];

% FFT
%----
Data = fft(Data);

% Filter
%-------
Data(FrequencySetZero_Index,:) = 0;

% IFFT
%-----
Data = ifft(Data);

% Undo zero padding
%------------------
Data_Filtered = Data(1:sampleLength,:);

end


function Result = nextpow2_one35(n)
% 
%Compute the min length for FFT according to AFNI's algorithm, By Xiao-Wei Song
%--------------------------------------------------------------------------
% 
% Original function name: rest_nextpow2_one35
% 
%	Copyright(c) 2007~2010
%	State Key Laboratory of Cognitive Neuroscience and Learning in Beijing Normal University
%	Written by Xiao-Wei Song 
%	http://resting-fmri.sourceforge.net
% 	<a href="Dawnwei.Song@gmail.com">Mail to Author</a>: Xiaowei Song
%	Version=1.0;
%	Release=20070903;

if length(n)>1
    n = cast(length(n),class(n));
end
if n<16
    Result =2^nextpow2(n);
    return;
end 

limit = nextpow2(n);           % n=134, limit=8
tbl   = [2^(limit-1):2^limit]; % tbl = 128, 129, ... , 256
tbl   = tbl(tbl>=n);           % tbl = 134, 135, ... , 256
for x=1:length(tbl)
    Result = tbl(x);
    [f,p]  = log2(Result);
    if ~isempty(f) && f == 0.5       % Copy from nextpow2.m
        return;
    end
    if mod(Result,3*5)==0        
        y    = Result /(3*5);
        [f,p]= log2(y);
        if ~isempty(f) && f == 0.5   % Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,3)==0        
        y     = Result /3;
        [f,p] = log2(y);
        if ~isempty(f) && f == 0.5   % Copy from nextpow2.m
            return;
        end
    end
    if mod(Result,5)==0        
        y     = Result /5;
        [f,p] =log2(y);
        if ~isempty(f) && f == 0.5   % Copy from nextpow2.m
            return;
        end
    end
end
Result = NaN; % Should not reach, except when n=1

% csfft_nextup35 in AFNI list 1~1024, 20070516, dawnsong
% 2
% 4
% 6
% 8
% 10
% 12
% 16
% 20
% 24
% 30
% 32
% 40
% 48
% 60
% 64
% 80
% 96
% 120
% 128
% 160
% 192
% 240
% 256
% 320
% 384
% 480
% 512
% 640
% 768
% 960
% 1024

end