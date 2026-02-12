function giniCoeff = perpl_ComputeGiniCoeff(vector)
if size(vector, 1) >  size(vector, 2)
    vector = vector';
end

vector = sort(abs(vector), 'ascend');

cum1 = cumsum(vector);
cum2 = cum1/max(cum1);
giniCoeff = 1-2*sum(cum2)./(length(cum2)+1);

end

