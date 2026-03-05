### RNAseq data processing ### 
# 1. Raw read preprocessing
FastQC quality check (before trimming)
fastqc -o ${output_folder} ${read_file}_1.fastq
fastqc -o ${output_folder} ${read_file}_2.fastq
Adapter and quality trimming with Trimmomatic (PE reads)
java -jar trimmomatic-0.38.jar PE \
  ${read_file}_1.fastq ${read_file}_2.fastq \
  ${read_file}_1_trim.fq ${read_file}_1_trim_unpaired.fq \
  ${read_file}_2_trim.fq ${read_file}_2_trim_unpaired.fq \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:150
FastQC quality check (after trimming)
fastqc -o ${output_folder} ${read_file}_1_trim.fq
fastqc -o ${output_folder} ${read_file}_2_trim.fq

---

# 2. Read mapping with HISAT2
Build reference index 
hisat2-build reference_genome.fasta reference_genome_index
Mapping paired-end reads to the reference genome
hisat2 -p 8 \
  --no-unal \
  --rna-strandness FR \
  -x reference_genome_index \
  -1 ${read_file}_1_trim.fq \
  -2 ${read_file}_2_trim.fq \
  -S ${read_file}.sam

---

# 3. Post-mapping filtering
Extract uniquely mapped reads
egrep '^@|NH:i:1' ${read_file}.sam > ${read_file}_uniq.sam
Filter reads with mismatches (allowing ≤2 mismatches)
egrep -w 'XM:i:0|XM:i:1|XM:i:2' ${read_file}_uniq.sam > ${read_file}_uniq_mis.sam
Generate final uniquely mapped and filtered SAM file
cat header.sam ${read_file}_uniq_mis_filtered.sam > ${read_file}_final.sam

---

# 4. Gene-level quantification with HTSeq
Sort SAM file by read name
samtools sort -n -O sam ${read_file}_final.sam > ${read_file}_final_sorted.sam
Count reads per gene
htseq-count \
  --order=name \
  --mode=union \
  --stranded=reverse \
  --idattr=gene_id \
  -t mRNA \
  -f sam \
  ${read_file}_final_sorted.sam \
  reference_annotation.gtf \
  > ${read_file}_counts.txt
