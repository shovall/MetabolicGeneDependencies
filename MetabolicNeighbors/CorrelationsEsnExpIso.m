function resExpEsnIso = CorrelationsEsnExpIso(resExpEsn,resDist,minDist,maxDist,boolUseTheGeneItself,numGenesRequired)

resExpEsnIso = struct;
resExpEsnIso.genesEsn = resExpEsn.genesEsn;
resExpEsnIso.genesExp = resExpEsn.genesExp;

resExpEsnIso.r = zeros(size(resExpEsn.r))+NaN;
resExpEsnIso.p = zeros(size(resExpEsn.r))+NaN;

for i=1:length(resExpEsnIso.genesEsn)
    gene = resExpEsnIso.genesEsn{i};
    [~,loc] = ismember(gene,resDist.genes);
    if(~isnan(loc))
        locsRelatedGenes = find(resDist.data{loc}.dist<=maxDist);
        
        [sorted,~] = sort(resDist.data{loc}.dist);
        n = min(length(sorted),numGenesRequired);
        thresh = max(sorted(1:n));
        
        thresh  = min(thresh,maxDist);
        locsRelatedGenes = find(resDist.data{loc}.dist<=thresh & resDist.data{loc}.dist>=minDist);
    else
        locsRelatedGenes = {};
    end
    
    if(boolUseTheGeneItself)
        relatedGenes = {gene};
    else
        relatedGenes = {};
    end
    if(~isempty(locsRelatedGenes))
        relatedGenesPartial = resDist.data{loc}.genes(locsRelatedGenes);
        relatedGenes = union(relatedGenesPartial,relatedGenes);
    end
    if(~isempty(relatedGenes))
        [exists,locs] = ismember(relatedGenes,resExpEsnIso.genesExp);
        locs = locs(exists);
        resExpEsnIso.r(i,locs) = resExpEsn.r(i,locs);
        resExpEsnIso.p(i,locs) = resExpEsn.p(i,locs);
    end
end

resExpEsnIso.pFDRBonferoni = calcFDR(resExpEsnIso.p,1);
resExpEsnIso.pFDR = calcFDR(resExpEsnIso.p,0);
end
