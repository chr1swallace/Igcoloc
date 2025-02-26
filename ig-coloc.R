#!/usr/bin/env Rscript

args=R.utils::cmdArgs(default=list(args=TRUE,
                                   infiles="/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz,/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
                                   iregions="data/iga.regions,data/igm.regions",
                                   names="igg,igm",
                                   outfile="data/coloc_igg_igm.csv",
                                   outplot="data/coloc_igg_igm.pdf" ))

print(args)
vargs <- lapply(args[-1], function(x) strsplit(x, ",")[[1]])
trimfiles = sub(".regions","_trimdata.rds",vargs$iregions)

library(gridExtra)
library(ggplot2); theme_set(theme_bw())
library(coloc)
library(data.table)
library(patchwork)

regions=lapply(vargs$iregions,fread)
names(regions)=vargs$names
for(i in 1:2) {
    regions[[i]]$ig=vargs$names[i]
    regions[[i]]=regions[[i]][!(chromosome==6 & base_pair_location > 25e+6 & base_pair_location < 34e+6)] # remove MHC
}

## find overlaps
library(magrittr)
library(GenomicRanges)
window=5e+5
gregions=lapply(regions, function(x) GRanges(seqnames=x$chromosome, ranges=IRanges(start=x$base_pair_location-window, end=x$base_pair_location+window)))
overlaps=GenomicRanges::findOverlaps(gregions[[1]],gregions[[2]])  %>%
    as.data.table()

## load data (slow)
data=lapply(trimfiles, readRDS) 
for(i in 1:2) {
    data[[i]][,snp:=paste(chromosome,base_pair_location,other_allele,effect_allele,sep="_")]
    data[[i]]=data[[i]][!duplicated(snp)]
}

## add sample size (also slow)
## ss=sub("meta.tsv.gz","per_snp_sample_size.tsv.gz",vargs$infiles)  %>%
##     lapply(., fread)
for(i in 1:2) {
    data[[i]]=merge(data[[i]],ss[[i]][,.(chromosome,base_pair_location,other_allele,effect_allele,sample_size)],
                    by=c("chromosome","base_pair_location","other_allele","effect_allele"))
    tt=table(data[[i]]$sample_size)
    data[[i]][,m3:=sample_size > as.numeric(names(tt[tt>3e+6]))]
    data[[i]][,sample_frac:=sample_size/max(sample_size, na.rm=TRUE)]
    data[[i]]=data[[i]][!duplicated(snp)]
}

## library(ggplot2)
## p1=ggplot(data[[1]], aes(x=sample_size)) + geom_histogram()
## p2=ggplot(data[[2]], aes(x=sample_size)) + geom_histogram()
## library(patchwork)
## p1 / p2

results=plots=vector("list",nrow(overlaps))
for(i in 1:nrow(overlaps)) {
    region1=regions[[1]][overlaps[i]$queryHits]
    region2=regions[[2]][overlaps[i]$subjectHits]
    st=min(region1$base_pair_location-window,region2$base_pair_location-window)
    en=max(region1$base_pair_location+window,region2$base_pair_location+window)
    x=data[[1]][chromosome==region1$chromosome & base_pair_location>=st & base_pair_location<=en]
    y=data[[2]][chromosome==region1$chromosome & base_pair_location>=st & base_pair_location<=en]
    m=merge(x,y, by=c("snp"), suffixes=c(".x",".y"))
    use=m[abs(sample_frac.x-sample_frac.y) < 0.1]$snp
    x=x[snp %in% use]
    y=y[snp %in% use]
    ## p1=ggplot(m, aes(x=sample_size.x)) + geom_histogram()
    ## p2=ggplot(m, aes(x=sample_size.y)) + geom_histogram()
    ## p3=ggplot(m, aes(x=sample_size.x,y=sample_size.y)) + geom_point()
    ## (p1 / p2) + p3
    ## with(m, table(m3.x,m3.y))

    result=coloc.abf(with(x[], list(beta=beta,
                          varbeta=standard_error^2,
                          snp=snp,
                          type="quant",
                          sdY=1)),
                     with(y[],list(beta=beta,
                          varbeta=standard_error^2,
                          snp=snp,
                          type="quant",
                          sdY=1)))

    m=merge(x[],y[],by=c("base_pair_location","snp"),all=TRUE,suffixes=c(".1",".2"))[,c("base_pair_location",paste0(c("p_value.","beta.","standard_error."),rep(1:2,each=3))),with=FALSE] 
    mm=melt(m, id.vars="base_pair_location", measure.vars=patterns("p_value","beta","standard_error"), variable.name="ig", value.name=c("p_value","beta","standard_error"))
    mm[,ig:=vargs$names[ig]]
    p.manh=ggplot(mm) +
        geom_point(aes(x=base_pair_location,y=-log10(p_value))) +
        facet_grid(ig ~ .)
    p.beta=ggplot(m) +
        geom_point(aes(x=beta.1,y=beta.2)) +
        labs(x=paste0("beta ",vargs$names[1]),y=paste0("beta ",vargs$names[2]))
    p.z=ggplot(m) +
        geom_point(aes(x=beta.1/standard_error.1,y=beta.2/standard_error.2)) +
        labs(x=paste0("z ",vargs$names[1]),y=paste0("z ",vargs$names[2]))
    plots[[i]]=p.manh + (p.beta / p.z) +
        plot_annotation(subtitle=paste0("chromosome ",region1$chromosome,":",region1$base_pair_location,",",region2$base_pair_location,"\t",paste(vargs$names,collapse=" vs ")),
                        title=paste0("pp.H4=",format.pval(result$summary["PP.H4.abf"]),"\t","pp.H3=",format.pval(result$summary["PP.H3.abf"])))
    results[[i]]=result$summary
}

results %<>% do.call("rbind",.)  %>%
    as.data.table()
results[,ig1:=vargs$names[1]]
results[,ig2:=vargs$names[2]]
results[,chromosome:=regions[[1]][overlaps$queryHits]$chromosome]
results[,location1:=regions[[1]][overlaps$queryHits]$base_pair_location]
results[,location2:=regions[[2]][overlaps$subjectHits]$base_pair_location]
results[,d:=abs(location2-location1)]
fwrite(results,args$outfile)

pdf(args$outplot, onefile = TRUE)
for (i in seq(length(plots))) {
  print(plots[[i]])  
}
dev.off()
