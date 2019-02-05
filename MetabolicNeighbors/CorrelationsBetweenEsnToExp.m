function [resExpEsn] = CorrelationsBetweenEsnToExp(depMapData,profileDataCCLE,genesSelected)

resExpEsn = struct;
resExpEsn.genesEsn = genesSelected;
resExpEsn.genesExp = profileDataCCLE.common;
r = zeros(length(genesSelected),length(resExpEsn.genesExp));
p = zeros(length(genesSelected),length(resExpEsn.genesExp))+1;

[exists,locsExp] = ismember(depMapData.celllines,profileDataCCLE.caseId);
locsEsn = find(exists);
locsExp = locsExp(locsEsn);

numGenesEsn = length(genesSelected);
numGesnesExp = length(resExpEsn.genesExp);

parfor i=1:numGenesEsn
    if(mod(i,50)==0)
        disp(i);
    end
    [~,locGene] = ismember(genesSelected{i},depMapData.genes);
    dataEsn = depMapData.data(locGene,locsEsn);
    for j=1:numGesnesExp
        dataExp = profileDataCCLE.data(j,locsExp);
        locs1 = find(~isnan(dataExp));
        locs2 = find(~isnan(dataEsn));
        locs = intersect(locs1,locs2);
        if(isempty(locs))
            r(i,j)  = NaN;
            p(i,j) = NaN;
            continue;
        end
        [rCurr,pCurr] = corr(dataEsn(locs)',dataExp(locs)','type','Spearman');
        r(i,j)  = rCurr;
        p(i,j) = pCurr;
    end
end
resExpEsn.r = r;
resExpEsn.p = p;
resExpEsn.pFDRBonferoni = calcFDR(resExpEsn.p,1);
resExpEsn.pFDR = calcFDR(resExpEsn.p,0);
end

