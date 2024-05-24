function ppc0 = perpl_ppc(crss)
% 'Unweighted PPC0'
% input [spikes freq/cmx values]
% Input is extracted complex values from TF based on spike times

crss = crss ./ abs(crss); % normalize the angles before averaging

%% PPC Compute pairwise phase consistency (PPC0 of Vinck et al, 2011)
dim = 1;
dof = sum(~isnan(crss),dim);
sinSum = abs(nansum(imag(crss),dim));
cosSum = nansum(real(crss),dim);
ppc0 = (cosSum.^2+sinSum.^2 - dof)./(dof.*(dof-1));  




