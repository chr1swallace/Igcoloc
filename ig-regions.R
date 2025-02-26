library(data.table)
library(coloc)
library(glue)
library(magrittr)

args=R.utils::cmdArgs(defaults=list(infile="/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
                                    outfile="data/regions_iga.tsv"))

## baseDir="/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline"
## f.iga=glue("{baseDir}/results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz")
## f.igg=glue("{baseDir}/results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz")
## f.igm=glue("{baseDir}/results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz")

## sample sizes
## iga: /rds/project/rdsHNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/per_snp_sample_size.tsv.gz
## igm: /rds/project/rdsHNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/per_snp_sample_size.tsv.gz
## igg: /rds/project/rdsHNdhZnUvWRk/analysis/pid/common_variant_analysis/serum_ig_pipeline/results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/per_snp_sample_size.tsv.gz
x=fread(args$infile)
ss=sub("meta.tsv.gz","per_snp_sample_size.tsv.gz",args$infile)  %>% fread()
x=merge(x,ss[,.(chromosome,base_pair_location,other_allele,effect_allele,sample_size)],
        by=c("chromosome","base_pair_location","other_allele","effect_allele"))
x[,sample_frac:=sample_size/max(sample_size, na.rm=TRUE)]
x=x[!(chromosome==6 & base_pair_location > 25e+6 & base_pair_location < 34e+6)] # remove MHC
x=x[!is.na(beta) & standard_error > 0] # remove SNPs with no effect size
x[,snp:=paste(chromosome,base_pair_location,other_allele,effect_allele,sep="_")]
x=x[!duplicated(snp)]

## find lead signals - lead snps +/- 1mb
y=x[p_value < 5e-8]
regions=vector("list",200) # oversize, can cut later
i=0
y$done=FALSE
window=5e+5
while(any(y[done==FALSE]$p_value < 5e-8)) {
    i=i+1
    w=which.min(y[done==FALSE]$p_value)
    message(i,"\t",y[done==FALSE]$p_value[w])
    regions[[i]]=y[done==FALSE][w,.(chromosome,base_pair_location)]
    y[chromosome==regions[[i]]$chromosome &
      (base_pair_location > regions[[i]]$base_pair_location - window) &
      (base_pair_location < regions[[i]]$base_pair_location + window), done:=TRUE]
}
regions  %<>% rbindlist()

fwrite(regions,file=args$outfile)

bigwindow=1.5e+6
x=x[chromosome %in% regions$chromosome]
x[,inregion:=FALSE]
for(i in 1:nrow(regions))
    x[ chromosome==regions$chromosome[i] &
       (base_pair_location > regions$base_pair_location[i] - bigwindow) &
       (base_pair_location < regions$base_pair_location[i] + bigwindow), inregion:=TRUE]
table(x$inregion)

x=x[inregion==TRUE]

saveRDS(x, file=sub(".regions","_trimdata.rds",args$outfile))
