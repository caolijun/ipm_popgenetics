####make folders in output folder of stacks
###mkdir genome populations stacks bam/sorted sam
cores=60
species=cs
samples=/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/1_indi_sort/                 ###reads data of individuals
out=/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/2_call_SNP/2_refer_based/bowtie  ###out put folders of stacks
genome=/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/2_call_SNP/2_refer_based/bowtie/cs.hic.genome.fasta ###reference genome
map=/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/2_call_SNP/2_refer_based/bowtie/cs.popmap.218.txt      ###popmap file for running populations

#######sorting individuals####

rawData=/ipm1/ljcao/cs.rad/1_raw_data

list="BJYQ HBYC LNXC HLHE HLMD SXJZ SXXY NXWZ SDTA SDLK YNKM SCAB"
for pop in $list
        do
reads1=${rawData}/${pop}.R1.fastq.gz
reads2=${rawData}/${pop}.R2.fastq.gz
barcode=/ipm1/ljcao/cs.rad/2_analysis_and_processing_of_data/1_indi_sort/barcode/${pop}.txt ###for individual sorting
process_radtags -1 ${reads1}/ -2 ${reads2} -i gzfastq -b ${barcode} -o ${samples} -E phred33 -c -q -r --inline_inline --renz_1 nlaIII --renz_2 aciI -D -s 20 -t 140 &> ${pop}.process.rad.log
done

########mapping by bowtie2#########
###build genome index by bowtie2###

bowtie2-build ${genome} ${out}/${species}.genome/${species}.genome

####using a traversal to mapping reads to the genome###

list="BJYQ01 BJYQ03 BJYQ05 BJYQ06 BJYQ07 BJYQ08 BJYQ09 BJYQ10 BJYQ11 BJYQ12 ...."
for i in $list
        do
bowtie2 -p 80  --ma 2 --mp 6,2 --np 1 --gbar 4 --rdg 5,5 --rfg 5,5 --local -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 --no-unal -x ${out}/${species}.genome/${species}.genome -1 $samples/"$i".1.fq.gz -2 $samples/"$i".2.fq.gz -U $samples/"$i".rem.1.fq.gz,$samples/"$i".rem.2.fq.gz -S $out/sam/"$i".sam
done

#####detect errors from log files###

grep -iE "\b(err|e:|warn|w:|fail|abort)" *.log > errors

###sam2bam and sort###

list="BJYQ01 BJYQ03 BJYQ05 BJYQ06 BJYQ07 BJYQ08 BJYQ09 BJYQ10 BJYQ11 BJYQ12 BJYQ27 BJYQ14 BJYQ15 BJYQ16 BJYQ19 BJYQ20 BJYQ22 BJYQ23 BJYQ24 BJYQ26 HBYC11 HBYC13 HBYC15 HBYC16 HBYC17 HBYC18 HBYC19 HBYC20 HBYC21 HBYC22 HBYC23 HBYC24 HBYC25 HBYC26 HBYC27 HBYC28 HBYC29 HBYC30 HBYC31 HBYC32 HLHE02 HLHE03 HLHE04 HLHE05 HLHE06 HLHE07 HLHE08 HLHE09 HLHE10 HLHE11 HLHE12 HLHE13 HLHE14 HLHE15 HLHE16 HLHE17 HLHE18 HLHE19 HLHE20 HLHE21 HLMD50 HLMD51 HLMD52 HLMD54 HLMD55 HLMD57 HLMD59 HLMD60 HLMD61 HLMD62 HLMD63 HLMD64 HLMD65 HLMD66 HLMD67 HLMD69 HLMD70 HLMD71 HLMD72 HLMD73 YNKM01 YNKM02 YNKM03 YNKM04 YNKM05 YNKM06 YNKM07 YNKM08 YNKM09 YNKM10 YNKM11 YNKM12 YNKM13 YNKM14 YNKM15 YNKM16 YNKM17 YNKM18 YNKM19 YNKM20"
for i in $list
 do
samtools view -@40 -bS ${out}/sam/"$i".sam > ${out}/bam/"$i".bam
samtools sort --threads 40 ${out}/bam/"$i".bam > ${out}/bam/sorted/"$i".bam
done

########call SNP#########

ref_map.pl -T 200 --samples ${out}/bam/sorted --popmap ${map} -o ${out}/stacks -X "populations: --plink --vcf --fstats -t 100 -r 0.8 -p 12 --min_mac 2 --max_obs_het 0.7"

######we can run populations sole after ref_map.pl finished.##

populations -P ${out}/stacks/ -M ${map} -O ${out}/populations --plink --vcf --fstats -t 100 -r 0.8 -p 12 --min_mac 2 --min_maf 0.05 --max_obs_het 0.7

