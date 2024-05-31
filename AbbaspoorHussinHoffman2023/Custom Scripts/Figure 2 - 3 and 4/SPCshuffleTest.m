function  surrogate_spk_LFP_Convol = SPCshuffleTest(decimateCSC)

rng('shuffle')
peusoSpikeTimes = randperm(length(decimateCSC), 100000);
psuedoSpikeLFP = zeros(2, length(decimateCSC));
psuedoSpikeLFP(1,:) = decimateCSC;
psuedoSpikeLFP(2, peusoSpikeTimes) = 1;

trials{1} = psuedoSpikeLFP;
times{1} = [0:0.001:(length(decimateCSC)/1000)-0.001];

surrogate_spk_LFP_Convol = spkspectrum(trials, times);

end