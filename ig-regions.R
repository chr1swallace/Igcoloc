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

x=fread(args$infile)

## find lead signals - lead snps +/- 1mb
y=x[p_value < 5e-8]
regions=vector("list",200) # oversize, can cut later
i=0
y$done=FALSE
window=1e+6
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
