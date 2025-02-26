#!/usr/bin/env rake
require 'rake'
require 'date'
require '/home/cew54/slurmer/Qsub.rb'

QFILE="rake-#{DateTime.now}.sh"
qfile_=open(QFILE,'w')
qrun=false
qstr=" -r " # q options
args = {#:job=>"impute",
  :time=>'4:00:00',
  :tasks=>1,
  :cpus=>1,
  :autorun=>true,
  :interactiverun=> ! ENV['INTERACT'].nil?,
  :excl=>" "}


baseDir="/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline"
inp={:iga=>"#{baseDir}/results/iga_meta/with_epic/with_liu/with_scepanovic/with_dennis/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
     :igg => "#{baseDir}/results/igg_meta/with_epic/with_dennis/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz",
      :igm => "#{baseDir}/results/igm_meta/with_epic/with_scepanovic/with_pietzner/without_gudjonsson/with_eldjarn/meta.tsv.gz" }

region_files=inp.keys.map do |k|
  "data/#{k}.regions"
end

## * regions

desc "create region files"
task :regions => region_files do
  qrun=true
  qstr+=" -t 00:59:00 -j igregions -g 3 "
end

region_files.each_with_index do |f,i|
  file region_files[i] do |t|
    qfile_.puts "Rscript ig-regions.R --args infile=#{inp[inp.keys[i]]} outfile=#{region_files[i]}"
  end
end

## * coloc

# awkward stuff to make all pairwise combinations
combs=(0..2).to_a.combination(2).to_a
keycombs=combs.map do |c| inp.keys.values_at(c[0], c[1]) end
outfiles=keycombs.map do |c| "data/coloc_#{c.join('_')}.csv" end
outplots=keycombs.map do |c| "data/coloc_#{c.join('_')}.pdf" end
infiles1=combs.map do |c| inp.values[c[0]] end
infiles2=combs.map do |c| inp.values[c[1]] end
rfiles1=combs.map do |c| region_files[c[0]] end
rfiles2=combs.map do |c| region_files[c[1]] end
outfiles.each_with_index do |f,i|
  file outfiles[i] do |t| 
    qfile_.puts "./ig-coloc.R --args infiles=\"#{infiles1[i]},#{infiles2[i]}\" iregions=\"#{rfiles1[0]},#{rfiles2[i]}\" names=\"#{keycombs[i].join(',')}\" outfile=#{outfiles[i]} outplot=#{outplots[i]} > data/logs/#{File.basename(outfiles[i],'.csv')}.log 2>&1"
  end
end 

desc "run coloc"
task :coloc => outfiles do
  qrun=true
  qstr+=" -j coloc -g 3 -t 00:59:00 "
end

## * disease coloc

disease_files=Dir.glob("/rds/project/who1000/rds-who1000-wgs10k/analysis/pid/common_variant_analysis/serum_ig_pipeline/resources/harmonised_gwas/*.gz").grep_v(/iga|igg|igm/)
disease_files.each do |f|
  inp.each do |k,v|
    outfile="data/coloc_#{File.basename(f,'.tsv.gz')}_#{k}.csv"
    outplot="data/coloc_#{File.basename(f,'.tsv.gz')}_#{k}.pdf"
    task :dcoloc => outfile
    file outfile do
      qfile_.puts "Rscript disease-coloc.R --args ig=#{k} rfile=data/#{k}.regions igfile=#{v} diseasefile=#{f} outfile=#{outfile} outplot=#{outplot} > data/logs/#{File.basename(outfile,'.csv')}.log 2>&1"
    end
  end
end

task :dcoloc => disease_files do
  qrun=true
  qstr+=" -j dcoloc -t 00:59:00 "
end

## * run    

at_exit do
  qfile_.close
  if qrun
    n=`cat #{QFILE} | wc -l`.chomp
    if n.to_i > 1000
      g=n.to_i / 1000 + 1 # +1 because 6/5=1
      puts "\nTo run #{n} jobs:\nqlines_asarray.rb -g #{g} #{qstr} #{QFILE}"
    else 
      puts "\nTo run #{n} jobs:\nqlines_asarray.rb #{qstr} #{QFILE}"
    end
  else
    File.unlink(QFILE)
  end
end

