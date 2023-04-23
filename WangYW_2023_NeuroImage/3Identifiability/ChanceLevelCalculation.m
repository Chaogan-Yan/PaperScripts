

nGroups=41;
nSamplesForEachGroup=3;
nSamples=nSamplesForEachGroup*nGroups;

CumulativeCorrectValue=0;
for iCorrectGroups=nGroups:-1:1
    if iCorrectGroups==nGroups
        CorrectValue = (iCorrectGroups/nGroups) * 1;
    else
        CorrectValue= (iCorrectGroups/nGroups) * (stirling2((nSamples-iCorrectGroups*nSamplesForEachGroup),nGroups-iCorrectGroups) - stirling2(nSamples-(iCorrectGroups+1)*nSamplesForEachGroup,nGroups-(iCorrectGroups+1))); 
    end

    CumulativeCorrectValue=CumulativeCorrectValue+CorrectValue;
end

ChanceLeval=double(CumulativeCorrectValue)/double(stirling2(nSamples,nGroups));

