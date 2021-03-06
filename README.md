# GenomeAssemblies


# 1. Prepare data

- Download data
From public repositories as SRA https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=search_obj

- Change Raw reads to format XXX_R1.fastq.gz


# 2. Trimm sequences

    a=0;for folder in $(ls -d */); do cd $folder;a=0;for i in $(ls *_R1.fastq.gz | sort -t'_' -k2); do a=$((a + 1)); b=$(echo $i | cut -d'.' -f1 | cut -d'_' -f1); c=$(echo $i | sed 's/R1/R2/g');echo $i" "$c ; bsub -q normal -L /bin/bash -J TRIMMO$a -u ivan.mateusgonzalez@epfl.ch  -N  "  java -jar /home/imateus2/software/Trimmomatic-0.36/trimmomatic-0.36.jar PE -phred33 $i $c $(echo $i | cut -d'_' -f1,2)_Trimm_Pair.fastq $(echo $i | cut -d'_' -f1,2)_Trimm_Unpair.fastq $(echo $c | cut -d'_' -f1,2)_Trimm_Pair.fastq $(echo $c | cut -d'_' -f1,2)_Trimm_Unpair.fastq LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:50";  done ; cd ..; done


# 3. Velvet de novo assembly  in cluster

- 1 sample

      # testing multiple Kmers from K=35 to K=96 each 4 K
      module add UHTS/Assembler/velvet/1.2.10; velveth out_temp_ALLK 17,96,4 -fastq.gz -shortPaired -separate DT7_R1.fastq.gz_Trimm_Pair.fastq DT7_R2.fastq.gz_Trimm_Pair.fastq

      # assembling
      module add UHTS/Assembler/velvet/1.2.10; velvetg out_temp_35 -clean yes -exp_cov auto -cov_cutoff auto -min_contig_lgth 200

- Evaluate parameters effect on genome assembly

      array=( 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200)
      for i in "${array[@]}"; do echo $i ; module add UHTS/Assembler/velvet/1.2.10; velvetg out_ALLK_25 -clean yes -exp_cov $i -cov_cutoff auto -min_contig_lgth 200; done | grep "Final graph" > OUT_25.txt

      # loop to do for all the kmers, and strains.
      # finding the good parameters Exp_cov

      Kmer=( 17 21 25 29 33 37 41 45 49 53 57 61 65 69 73 77 81 85 89 93)
      array=( 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200)

      a=0;for folder in $(ls -d *_ALLK_*/); do echo $(echo $folder | sed 's/\///') ; done

      for k in "${Kmer[@]}"; do echo out_temp_ALLK_$(echo $k) ;for i in "${array[@]}"; do echo $i ; module add UHTS/Assembler/velvet/1.2.10; velvetg out_temp_ALLK_$(echo $k) -clean yes -exp_cov $i -cov_cutoff auto -min_contig_lgth 200; done ; done| grep "Final graph" > OUT_1992.txt

- Choose parameter per sample with script Estimation_Velvet_assemblies_parameters_v1.R

-  building brujin graphh kmers 

        a=0;for folder in $(ls -d Th*/); do cd $folder;a=0;for l in $(ls *_R1.*Pair.fastq | sort -t'_' -k2); do a=$((a + 1)); b=$(echo $l | cut -d'.' -f1 | cut -d'_' -f1); m=$(echo $l | sed 's/R1/R2/g');echo $l" "$m ; bsub -q normal -L /bin/bash -J velveth$a -u ivan.mateusgonzalez@epfl.ch  -N -R "rusage[mem=16000]" -M 16000000 " module add UHTS/Assembler/velvet/1.2.10; velveth $(echo $l | cut -d'_' -f1)_ALLK 17,96,4 -fastq.gz -shortPaired -separate $l $m " ; done ; cd .. ; done 


- Assembling

        a=0;for folder in $(ls -d *_ALLK_*/); do echo $(echo $folder | sed 's/\///') ;array=(10 20 30 40 50 60 70 80 90 100 110 120 130 140 150 160 170 180 190 200); for i in "${array[@]}"; do echo $i ;  bsub -q normal -L /bin/bash -J velvetg_$folder -u ivan.mateusgonzalez@epfl.ch  -N -R "rusage[mem=16000]" -M 16000000  "module add UHTS/Assembler/velvet/1.2.10 ; velvetg $folder -clean yes -exp_cov $i -cov_cutoff auto -min_contig_lgth 200 > OUT_$(echo $folder | sed 's/\///')_$i.txt "; done  ; done


- get values
            
         for i in $(ls OUT*.txt); do cat $i  ; done | grep "Final graph" | sed 's/Final graph has //' | sed 's/nodes and n50 of //'  | sed 's/, max//'  | sed 's/, total//'  | sed 's/, using//' | sed 's/ reads//' | sed 's/\// /'

- Circularization

      module add UHTS/Analysis/amos/3.1.0;
      toAmos -s contigs.fa -o circularized.afg
      minimus2 circularized

      bsub -q normal -L /bin/bash -J min -u ivan.mateusgonzalez@epfl.ch  -N  " module add UHTS/Analysis/amos/3.1.0; toAmos -s contigs.fa -o circularized.afg  ; minimus2 circularized "


