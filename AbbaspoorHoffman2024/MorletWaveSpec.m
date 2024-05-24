%%% Modified from [wavespec] = bz_WaveSpec(lfp) calculates the
%%% Saman Abbaspoor - 03/12/2024 - Hoffman Lab - Vanderbilt

function [wavespec] = MorletWaveSpec(lfp,varargin)
%[wavespec] = bz_WaveSpec(lfp) calculates the
%wavelet transform of a signal with nfreqs frequencies in the range frange
%[fmin fmax]. Spacing between frequencies can be 'lin' or 'log'.
%Time-frequency resolution is defined by ncyc, the number of cycles in each
%wavelet. Uses Morlet (Gabor) wavelet.
%
%
%INPUT
%    lfp            a buzcode structure with fields lfp.data,
%                                                   lfp.timestamps
%                                                   lfp.samplingRate
%                   -lfp can also be a [t x 1] timeseries signal. in which
%                   case you need to input 'samplingRate'
%    <options>      optional list of property-value pairs (see table below)
%
%    =========================================================================
%     Properties    Values
%    -------------------------------------------------------------------------
%       'frange'    [low frequency, high frequency]     (default: [1 128])
%       'nfreqs'    number of frequencies               (default: 100
%       'roundfreqs' round freqs to unique integer vals (default: false
%                       *Note this may decrease number of freqs
%       'ncyc'      number of cycles in your wavelet    (default: 5)
%       'fvector'   predefined vector of frequencies
%       'space'     'log' or 'lin'  spacing of f's      (default: 'log')
%       'samplingRate' (only if input is not a buzcode structure)
%       'intervals'  restrict your spectrogram to timestamps in specific
%                   intervals
%       'showprogress' true/false (default:false)
%    =========================================================================
%
%OUTPUT
%   wavespec            buzcode-style structure
%       .data           [t x nfreqs] your spectrogram
%       .timestamps     [t x 1] timestamps
%       .freqs          frequencies of each column
%       .samplingRate   (Hz)
%       .channels       Channel indices of channels filtered... taken from lfp input
%       .filterparms    a structure that holds the parameters used for
%                       filtering, for future reference
%

%Last Updated: 10/9/15
%DLevenstein
%Modified by Antonio FR, 7/18/18
%Modified 4/30/2022

%% Parse the inputs

%Parameters
parms = inputParser;
addParameter(parms,'frange',[1 128],@isnumeric);
addParameter(parms,'nfreqs',100,@isnumeric);
addParameter(parms,'ncyc',5,@isnumeric);
addParameter(parms,'space','lin');
addParameter(parms,'samplingRate',[]);
addParameter(parms,'showprogress',false,@islogical);
addParameter(parms,'roundfreqs',false,@islogical);
addParameter(parms,'fvector',[]);

parse(parms,varargin{:})
frange = parms.Results.frange;
nfreqs = parms.Results.nfreqs;
ncyc = parms.Results.ncyc;
space = parms.Results.space;
samplingRate = parms.Results.samplingRate;
showprogress = parms.Results.showprogress;
roundfreqs = parms.Results.roundfreqs;
fvector = parms.Results.fvector;


%lfp input
if isstruct(lfp)
    samplingRate = lfp.samplingRate;
elseif isempty(lfp)
    wavespec = lfp;
    return
elseif isnumeric(lfp)
    data_temp = lfp;
    clear lfp
    lfp.data = data_temp;
    lfp.timestamps = [1:length(lfp.data)]'./samplingRate;
end

si = 1./samplingRate;

%%
if ~isa(lfp.data,'single') || ~isa(lfp.data,'double')
    lfp.data = single(lfp.data);
end

%Frequencies
if ~isempty(fvector)
    freqs = fvector;
else
    fmin = frange(1);
    fmax = frange(2);
    if strcmp(space,'log')
        assert(fmin~=0,'Log-spaced frequencies cannot have min of 0')
        freqs = logspace(log10(fmin),log10(fmax),nfreqs);
    elseif strcmp(space,'lin')
        freqs = linspace(fmin,fmax,nfreqs);
    else
        display('Frequency spacing must be "lin" or "log".')
    end
end

if roundfreqs
    freqs = unique(round(freqs));
end

%Filter with wavelets
nfreqs = size(freqs,2);
nchan = size(lfp.data,2);
ntime = ceil(size(lfp.data,1));
wavespec.data = nan(ntime,nfreqs,nchan);
wavespec.timestamps = lfp.timestamps;
for cidx = 1:nchan
    for f_i = 1:nfreqs
        if showprogress
            bz_Counter(f_i,nfreqs,'Wavelet Frequency')
        end
        wavelet = MorletWavelet(freqs(f_i),ncyc,si);
        wavespec.data(:,f_i,cidx) = ...
            FConv(wavelet',lfp.data(:,cidx));
    end
end

%% Output
wavespec.freqs = freqs;
wavespec.nfreqs = nfreqs;
wavespec.samplingRate = samplingRate;
wavespec.filterparms.ncyc = ncyc;
wavespec.filterparms.nfreqs = nfreqs;
wavespec.filterparms.frange = frange;
wavespec.filterparms.space = space;

clear lfp

end



function [wavelet, t] = MorletWavelet( f, numcyc, si )

% parameters...
s = numcyc/(2*pi*f);    %SD of the gaussian
tbound = (4*s);   %time bounds - at least 4SD on each side, 0 in center
tbound = si*ceil(tbound/si);
t = -tbound:si:tbound;    % time

% and together they make a wavelet
sinusoid = exp(2*pi*1i*f.*t);
gauss = exp(-(t.^2)./(2*s^2));
A = 1;
wavelet = A * sinusoid .* gauss;

%normalize to unit energy
wavelet = wavelet./norm(wavelet);

end


function [ convolution_result_fft ] = FConv( kernel, signal, varargin )
%FConv(kernel,signal) convolves a signal with a kernal via the fourier
%transform for speedy delivery.
%Note: as is, kernel is centered.
%   Adapted from Cohen Chapter 13
%Last Updated: 7/30/15
%DLevenstein

% Kernel and signal must be same orientation


% FFT parameters
n_kernel            = length(kernel);
n_signal               = length(signal);
n_convolution        = n_kernel + n_signal-1;
half_of_kernel_size = (n_kernel-1)/2;


% FFT of wavelet and EEG data
fft_kernel = fft(kernel,n_convolution);
fft_signal = fft(signal,n_convolution);

convolution_result_fft = ifft((fft_kernel.*fft_signal),n_convolution);

% cut off edges
convolution_result_fft = convolution_result_fft(half_of_kernel_size+1:end-half_of_kernel_size);


end