## Pangenome graph building with pggb 

## Resources:
- [Github page](https://github.com/pangenome/pggb)
- [pggb nextflow](https://nf-co.re/pangenome/1.0.0/)

## Run pggb using nextflow nf-core
I could not get pggb to work when I downloaded it using conda, singluarity, manually, etc. Nextflow was the only thing to get the graph building step to work for me. I must have some conflict with some dependencies. However, pggb downloaded with conda worked for the initial partitioning steps, but the graph building steps (wfmash, seqwish, smoothxg, and gfaffix) did not work for me unless I used nf-core.
```bash
# Install with conda
module load python/3.10.9-fasrc01
conda activate env_nf
```

## 1. First we need to partition sequences into communities. I followed [this](https://gtpb.github.io/CPANG22/pages/Day2a_Homo_sapiens_pangenome_graphs) tutorial first to partition my sequences. 

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
cd ${PATH_SAMPLE_FA_GZ}

DIR_PARTITIONING="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb/partitioning"

mkdir -p parts

( seq 34 ) | while read i; do
  awk '$6 ~ "scaffold_'$i'$"' $(ls $DIR_PARTITIONING/*.vs.ref.paf | \
    grep -v unaligned | sort; echo $DIR_PARTITIONING/rescues.paf) | cut -f 1 | sort | uniq \
    > parts/scaffold_$i.contigs;
done

( seq 34 ) | while read i; do
  echo scaffold_$i
  samtools faidx HemMar_prefixed_scaffolds_corrected.fa.gz HemMar#1#scaffold_$i > parts/scaffold_$i.pan.fa

  ls *.p_ctg.fa.gz | while read FASTA; do
    NAME=$(basename $FASTA *.p_ctg.fa.gz);
    echo scaffold_$i $NAME

    samtools faidx $FASTA $( comm -12 <(cut -f 1 $FASTA.fai | sort) <(sort parts/scaffold_$i.contigs) ) >> parts/scaffold_$i.pan.fa
  done

  bgzip -@ 8 parts/scaffold_$i.pan.fa && samtools faidx parts/scaffold_$i.pan.fa.gz
done
```

## 2. Run pggb using nf-core

```bash
conda activate env_nf

in_dir="/n/holyscratch01/edwards_lab/Users/kelsielopez/hap_assemblies/prefixed/parts"
outdir="/n/holyscratch01/edwards_lab/Users/kelsielopez/pggb"

cd /n/holyscratch01/edwards_lab/Users/kelsielopez/pggb

ls $in_dir/*.pan.fa.gz | while read FASTA; do
  NAME=${FASTA%.pan.fa.gz};
  echo $NAME
  mkdir $NAME
  cd $NAME
    nextflow run nf-core/pangenome -r dev --input $FASTA --n_haplotypes 12 --wfmash_map_pct_id 85 --seqwish_min_match_length 79 --outdir $NAME -profile singularity
  cd ../
  done
```

## 3. Merge vcfs, clean, etc. 

## 4. Compute variant table

Without repeat information or genomic overlaps:

With repeat information and genomic overlaps:
