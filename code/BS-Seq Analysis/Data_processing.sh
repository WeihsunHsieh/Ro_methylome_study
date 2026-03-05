### Data processing of BS-Seq data ###
# 1.	Fastqc v0.11.8
  CMD: 
fastqc -o ${output_folder} ${read_file_name}_1.fq
fastqc -o ${output_folder} ${read_file_name}_2.fq

---

# 2.	Trim v0.38
	CMD:
java -jar trimmomatic-0.38.jar SE ${read_file_name}_1.fq ${read_file_name}_1.trim.fq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:150
java -jar trimmomatic-0.38.jar SE ${read_file_name}_2.fq ${read_file_name}_2.trim.fq \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:150

---

# 3.	bs-seeker2 build index: bs-seeker2 v2.1.8; bowtie2 v2.3.2
	CMD: bs_seeker2-build.py -p ${bowtie2_path} -d ${ reference_genome_Index_folder } \
--aligner=bowtie2 -f ${reference_genome}

---

# 4.	Mapping
	CMD: 
bs_seeker2-align.py --bt2-p 1 -m 0.05 -p ${bowtie2_path} -d ${reference_genome_Index_folder} \
-i ${read_file_name}_1.trim.fq -e 150 --aligner=bowtie2 -o ${read_file_name}_1.trim.bsker1 \
-f bs_seeker1 -g ${reference_genome}
bs_seeker2-align.py --bt2-p 1 -m 0.05 -p ${bowtie2_path} -d ${reference_genome_Index_folder} \
-i ${read_file_name}_2.trim.fq -e 150 --aligner=bowtie2 -o ${read_file_name}_2.trim.bsker1 \
-f bs_seeker1 -g ${reference_genome}
cat ${read_file_name}_1.trim.bsker1 ${read_file_name}_2.trim.bsker1 | sort -k1,1 > \
${read_file_name}.trim.bsker1

---

# 5.	Convert scaffold names in the mapping results
	Script: convert_BSseeker_scaffold_names.pl
	CMD: perl convert_BSseeker_scaffold_names.pl \
${species}_conversion_scaffold_names_annotated.txt ${read_file_name}.trim.bsker1 > \ ${read_file_name}.trim.bsker1.converted
	Format: ${species}_conversion_scaffold_names_annotated.txt (delimiter: space)
fasta abbrev length type
9999 Lambda 48502 Lambda
9996 Chloroplast 134546 Chloroplast
1 Chr1 45524768 Chromosome

---

# 6.	Remove the clonal reads by checking the unique mapped read separately for different libraries
	CMD: sort -k1,1n ${read_file_name}.trim.bsker1.converted | awk '!x[$0]++' > \ ${read_file_name}.trim.bsker1.converted.declonal

---

# 7.	Filtering non-conversion reads
	CMD: awk '{if($8==0) print $0}' ${read_file_name}.trim.bsker1.converted.declonal > \ ${read_file_name}.trim.bsker1.converted.declonal.filtered

---

# 8.	Categorizing the reads
	CMD: 
awk '{if($9==${category}) print $0}' ${read_file_name}.trim.bsker1.converted.declonal.filtered > \ ${read_file_name}.trim.bsker1.converted.declonal.filtered.${category}
chr_nums=(${chromosome//"-"/" "})
chr_start=${chr_nums[0]}
chr_end=${chr_nums[1]}
for i in `seq ${chr_start} ${chr_end}`
do
    awk '{if($9=="Chromosome"&&$1=="'${i}'") print $0}' \
${read_file_name}.trim.bsker1.converted.declonal.filtered > \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${i}
done
cat ${read_file_name}.trim.bsker1.converted.declonal.filtered.chr* > \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}

---

# 9.	Generate perC file
	Script: count_Methylation_perC_converted.pl
	CMD:
tmp=$(ls -1 ${read_file_name}.trim.bsker1.converted.declonal.filtered.*)
for i in ${tmp}
do
	perl count_Methylation_perC_converted.pl ${i}
done
cat ${read_file_name}.trim.bsker1.converted.declonal.filtered.chr*.perC > \ ${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}.perC

---

# 10.	Binomial processing
	Script: Calculate_Binomial.pl
	CMD: perl Calculate_Binomial.pl \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}.perC 0.005 > \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}.perC.binomial

---

# 11.	Filtering at least two reads per C site
	CMD: 
awk '{if ($5>=2||($5==0&&$6>=2)) print $0}' \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}.perC.binomial | \
sort -k1,1n > \
${read_file_name}.trim.bsker1.converted.declonal.filtered.chr${chromosome}.perC.binomial.filtered
