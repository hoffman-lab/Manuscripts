function bic = bicoher_montecarlo(signal_length, Freq, n_permutes,alpha)

h = waitbar(0,'Monte Carlo Significance');

bic = NaN(155, 155, n_permutes);

for ii=1:n_permutes
    
    RedNoise = rednoise(signal_length, 1);
    [tmp_bic, waxis] = bicoherence(RedNoise, 'nsamp', 1024);
    
    freq = waxis.*1000;
    [c, f1] = min(abs(freq - Freq(1)));
    [c, f2] = min(abs(freq - Freq(2)));
    tmp_bic = tmp_bic(f1:f2, f1:f2);
    bic(:, :, ii) = tmp_bic;
    waitbar(ii/n_permutes,h);
    
end
close(h);
% ptile = prctile(bic,100*(1-alpha),3);

end




