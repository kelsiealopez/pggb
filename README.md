## Pangenome graph building with pggb 

## Table of Contents
- [Github page](addlinkhere)
- [pggb nextflow](https://nf-co.re/pangenome/1.0.0/)

## Installation
Steps to install and set up the project. Include code snippets for installation commands.


I ran pggb using nextflow
```bash
# Install with conda
module load python/3.10.9-fasrc01
conda activate env_nf
```

First we need to partition sequences into communities. I followed [this](https://gtpb.github.io/CPANG22/pages/Day2a_Homo_sapiens_pangenome_graphs) tutorial first to partition my sequences. 

Partition contigs by chromosome by mapping each assembly against the scaffolded references:
```bash
# Using -s 5k -p 85 -N -m -t 8
# i used -p 85 because my outgroup is 14% divergent from Hemitriccus

PATH_REFERENCE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/HemMar_prefixed_scaffolds_final.fa.gz"
PATH_SAMPLE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

mkdir -p /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning

cd ${PATH_SAMPLE_FA_GZ}

ls *.p_ctg.fa.gz | while read FASTA; do
  NAME=$(basename $FASTA .p_ctg.fa.gz);
  echo $NAME

  PATH_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.paf"
  wfmash $PATH_REFERENCE_FA_GZ ${PATH_SAMPLE_FA_GZ}/$FASTA -s 5k -p 85 -N -m -t 8 > $PATH_PAF
done
```
Collect unmapped contigs and remap them in split mode:
```bash
# Using -s 5k -p 85 -N -m -t 8
PATH_REFERENCE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/HemMar_prefixed_scaffolds_final.fa.gz"
PATH_SAMPLE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

cd ${PATH_SAMPLE_FA_GZ}

ls *.p_ctg.fa.gz | while read FASTA; do
  NAME=$(basename $FASTA .p_ctg.fa.gz);
  echo $NAME

  PATH_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.paf"
  comm -23 <(cut -f 1 $FASTA.fai | sort) <(cut -f 1 $PATH_PAF | sort) > $NAME.unaligned.txt

  wc -l $NAME.unaligned.txt
  if [[ $(wc -l $NAME.unaligned.txt | cut -f 1 -d\ ) != 0 ]];
  then 
    samtools faidx $FASTA $(tr '\n' ' ' < $NAME.unaligned.txt) > $NAME.unaligned.fa
    samtools faidx $NAME.unaligned.fa

    PATH_NO_SPLIT_PAF="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning/$NAME.vs.ref.no_split.paf"
    wfmash $PATH_REFERENCE_FA_GZ ${PATH_SAMPLE_FA_GZ}/$NAME.unaligned.fa -s 5k -p 85 -m -t 8 > $PATH_NO_SPLIT_PAF
  fi
done
```
Collect our best mapping for each of our attempted split rescues:
```bash
cd /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning

ls *.vs.ref.no_split.paf | while read PAF; do
  cat $PAF | awk -v OFS='\t' '{ print $1,$11,$0 }' | sort -n -r -k 1,2 | \
    awk -v OFS='\t' '$1 != last { print($3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15); last = $1; }'
done > rescues.paf

```
Subset by chromosome, including the references
```bash

PATH_SAMPLE_FA_GZ="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed"

cd ${PATH_SAMPLE_FA_GZ}

DIR_PARTITIONING="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning"

mkdir -p parts

( seq 35 ) | while read i; do
  awk '$6 ~ "scaffold_'$i'$"' $(ls $DIR_PARTITIONING/*.vs.ref.paf | \
    grep -v unaligned | sort; echo $DIR_PARTITIONING/rescues.paf) | cut -f 1 | sort | uniq \
    > parts/scaffold_$i.contigs;
done

( seq 35 ) | while read i; do
  echo scaffold_$i
  samtools faidx HemMar_prefixed_scaffolds_final.fa.gz HemMar#1#scaffold_$i > parts/scaffold_$i.pan.fa

  ls *.p_ctg.fa.gz | while read FASTA; do
    NAME=$(basename $FASTA *.p_ctg.fa.gz);
    echo scaffold_$i $NAME

    samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort parts/scaffold_$i.contigs) ) >> parts/scaffold_$i.pan.fa
  done

  bgzip -@ 8 parts/scaffold_$i.pan.fa && samtools faidx parts/scaffold_$i.pan.fa.gz
done

``
