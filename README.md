
### Per sample analysis
```
NSLOTS=16
sample=sample_name
fq1=/path/to/file_R1.fq.gz
fq2=/path/to/file_R2.fq.gz
trimmed_fq1=/path/to/trimmed_R1.fq.gz
trimmed_fq2=/path/to/trimmed_R2.fq.gz
gtf=/path/to/gene_model.gtf
star_index=/path/to/star_index


# fastp, version 0.20.1
fastp \
-i ${fq1} \
-o ${trimmed_fq1} \
-I ${fq2} \
-O ${trimmed_fq2} \
-l 20 \
-3 -W 4 -M 20 \
-t 1 -T 1 \
-x \
--compression 1 \
--thread ${NSLOTS} \
--html ./fastp_html/${sample}.html


# mapping by STAR, version 2.7.11b
STAR --runThreadN ${NSLOTS} \
--genomeDir ${star_index} \
--readFilesIn ${trimmed_fq1} ${trimmed_fq2} \
--outFileNamePrefix ./${sample}. \
--readFilesCommand zcat \
--sjdbGTFfile ${gtf} \
--outSAMtype BAM Unsorted \
--outSAMunmapped Within KeepPairs \
--outSAMattributes NH AS nM NM jM ch cN \
--chimMultimapNmax 20 \
--chimSegmentMin 20 \
--chimJunctionOverhangMin 20 \
--chimMultimapScoreRange 3 \
--chimScoreJunctionNonGTAG -4 \
--chimNonchimScoreDropMin 10 \
--alignInsertionFlush Right \
--twopassMode Basic \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin 8 \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000

mv ${sample}.Aligned.out.bam ${sample}.bam


# markdup, samtools version 1.19.2
samtools fixmate -@ ${NSLOTS} -m -r ${sample}.bam ${sample}.fixmate.bam
samtools sort -@ ${NSLOTS} ${sample}.fixmate.bam -o ${sample}.sorted.bam
samtools markdup -@ ${NSLOTS} ${sample}.sorted.bam ${sample}.markdup.bam
samtools sort -@ ${NSLOTS} -n ${sample}.markdup.bam -o ${sample}.markdup.namesorted.bam


# detect chimera, pysam version 0.22.0
python detect_chimera_spilicing.py ${sample}.markdup.namesorted.bam
```
  
  
### Merge multiple samples
```
# merge samples
python merge_chimera_spilicing.py

# make metadata
python summarize_metadata.py
```
