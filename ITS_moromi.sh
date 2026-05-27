[TOC]

# 0. Set environment variables
    db=/opt/database/amplicon
    wd=~/project/amplicon/ITS_moromi
    sw=/opt/bin
    cd ${wd}

    conda activate amplicon
    mkdir -p seq result temp
    
# 1. ITS pipeline

## 1.1. Quality Control
    mkdir -p result/01_raw_qc
    fastqc seq/*.fq -o result/01_raw_qc -t 10
    
    mkdir -p result/02_raw_multiqc
    multiqc -f -d result/01_raw_qc -o result/02_raw_multiqc
    
    # Filter
    mkdir -p result/03_raw_fastp
    for i in `cat result/metadata.tsv |tail -n+2 |cut -f1`;do
      fastp -i seq/${i}_1.fq -I seq/${i}_2.fq \
        -o result/03_raw_fastp/${i}_1_fastp.fq -O result/03_raw_fastp/${i}_2_fastp.fq -w 10\
        -z 4 -q 20 -u 30 -n 10 -h result/03_raw_fastp/${i}_fastp.html -j result/03_raw_fastp/${i}_fastp.json 
    done
    
    for i in `cat result/metadata.tsv |tail -n+2 |cut -f1 `;do 
      cat result/03_raw_fastp/${i}_fastp.json |sed -n '6,25p' \
        |sed 's/\t//g'|sed 's/^.*://g' |tr '\n' '\t' \
        |awk 'BEGIN{OFS=FS="\t"}{print $1,$12,$16,$17,$18,$20}' \
        |sed 's/,\t/\t/g' |awk 'BEGIN{OFS=FS="\t"}{Raw=$1/2;Clean=$2/2;Q20=$3*100;Q30=$4*100;GC=$6*100;Effective=$2/$1*100; print Raw,Clean,Q20,Q30,$5,GC,Effective }' \
        |sed "s/^/${i}\t/" 
    done |sed '1i Sample ID\tRaw Reads\tClean Reads\tQ20(%)\tQ30(%)\tAvgLen(bp)\tGC(%)\tEffective(%)' \
    > result/03_raw_fastp/summary_qc.txt



## 1.2. Merge paired reads and label samples
    
    for i in `tail -n+2 result/metadata.tsv | cut -f 1`;do
      vsearch --fastq_mergepairs result/03_raw_fastp/${i}_1_fastp.fq --reverse result/03_raw_fastp/${i}_2_fastp.fq \
      --fastqout temp/${i}.merged.fq --relabel ${i}.
    done &
    
  
    # Merge all samples into a single file
    cat temp/*.merged.fq > temp/all.fq
  

## 1.3. Cut primers and quality filter

   vsearch --fastx_filter temp/all.fq \
      --fastq_stripleft 0 --fastq_stripright 0 \
      --fastq_maxee_rate 0.01 \
      --fastaout temp/filtered.fa


## 1.4 Dereplication
   vsearch --derep_fulllength temp/filtered.fa \
      --output temp/uniques.fa --relabel Uni --minuniquesize 40 --sizeout
    

## 1.5 Denoise ASV
    usearch -unoise3 temp/uniques.fa \
      -zotus temp/zotus.fa
    
    # Change the sequence name from Zotu to ASV for easier identification
    sed 's/Zotu/ASV_/g' temp/zotus.fa > temp/otus.fa

    mkdir -p result/04_raw_output
    cp -f temp/otus.fa result/04_raw_output/otus.fa


## 1.6. Creat feature table
    vsearch --usearch_global temp/filtered.fa --db result/04_raw_output/otus.fa \
    	--otutabout result/04_raw_output/otutab.txt --id 0.97 --threads 20
    

## 1.7 Species annotation – excluding plastids
    vsearch --sintax result/04_raw_output/otus.fa --db ${db}/usearch/unite-8.2.fasta \
      --tabbedout result/04_raw_output/otus.sintax --sintax_cutoff 0.8

    Rscript ~/scripts/otutab_filter_nonFungi.R \
      --input result/4_raw_output/otutab.txt \
      --taxonomy result/4_raw_output/otus.sintax \
      --output result/otutab.txt\
      --stat result/4_raw_output/otutab_nonFungi.stat \
      --discard result/4_raw_output/otus.sintax.discard

    # Sequences corresponding to the filter feature table
    cut -f 1 result/otutab.txt | tail -n+2 > result/otutab.id
    usearch -fastx_getseqs result/04_raw_output/otus.fa \
        -labels result/otutab.id -fastaout result/otus.fa
    
    # Filter feature table corresponding to sequence annotations
    awk 'NR==FNR{a[$1]=$0}NR>FNR{print a[$1]}'\
        result/04_raw_output/otus.sintax result/otutab.id \
        > result/otus.sintax
    
    # Fill in the last column
    sed -i 's/\t$/\td:Unassigned/' result/otus.sintax

    usearch -otutab_stats result/otutab.txt \
      -output result/otutab.stat
    
    cat result/otutab.stat
    # View the detailed data volume of the sample for resampling

## 1.8 Normlize by subsample
    mkdir -p result/05_alpha
    Rscript ~/scripts/otutab_rare.R --input result/otutab.txt \
      --depth 72699 --seed 123 \
      --normalize result/otutab_rare.txt \
      --output result/05_alpha/vegan.txt
    
    usearch -otutab_stats result/otutab_rare.txt \
      -output result/otutab_rare.stat
    

## 1.9 Alpha diversity

### 1.9.1. Calculate alpha diversity index
    usearch -alpha_div result/otutab_rare.txt \
      -output result/05_alpha/alpha.txt

### 1.9.2. Rarefaction
    time usearch -alpha_div_rare result/otutab_rare.txt \
      -output result/05_alpha/alpha_rare.txt -method without_replacement

### 1.9.3. Select high-abundance bacteria from each group for comparison
    Rscript ~/scripts/otu_mean.R --input result/otutab.txt \
      --design result/metadata.tsv \
      --group Stage --thre 0 \
      --output result/otutab_mean.txt

    # using an average abundance frequency of more than one in a thousand (0.1%)
    awk 'BEGIN{OFS=FS="\t"}{if(FNR==1) {for(i=3;i<=NF;i++) a[i]=$i;} \
        else {for(i=3;i<=NF;i++) if($i>0.1) print $1, a[i];}}' \
        result/otutab_mean.txt > result/05_alpha/otu_group_exist2.txt
  

## 1.10. Beta diversity
    mkdir -p result/06_beta/
    usearch -cluster_agg result/otus.fa -treeout result/otus.tree
    usearch -beta_div result/otutab_rare.txt -tree result/otus.tree \
      -filename_prefix result/06_beta/ 


## 1.11. Taxonomy summary
    cut -f 1,4 result/otus.sintax \
      |sed 's/\td/\tk/;s/:/__/g;s/,/;/g;s/"//g;s/\/Chloroplast//' \
      > result/taxonomy2.txt
    head -n3 result/taxonomy2.txt

   awk 'BEGIN{OFS=FS="\t"}{delete a; a["k"]="Unassigned";a["p"]="Unassigned";a["c"]="Unassigned";a["o"]="Unassigned";a["f"]="Unassigned";a["g"]="Unassigned";a["s"]="Unassigned";\
      split($2,x,";");for(i in x){split(x[i],b,"__");a[b[1]]=b[2];} \
      print $1,a["k"],a["p"],a["c"],a["o"],a["f"],a["g"],a["s"];}' \
      result/taxonomy2.txt > temp/otus.tax
    sed 's/;/\t/g;s/.__//g;' temp/otus.tax|cut -f 1-8 | \
      sed '1 s/^/OTUID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\n/' \
      > result/taxonomy.txt
    head -n3 result/taxonomy.txt

    mkdir -p result/07_tax
    for i in p c o f g s;do
      usearch -sintax_summary result/otus.sintax \
      -otutabin result/otutab_rare.txt -rank ${i} \
      -output result/07_tax/sum_${i}_all.txt
      
      cat result/07_tax/sum_${i}_all.txt |awk 'BEGIN{OFS=FS="\t"} {$NF="";print $0}' |sed 's/\t$//' > result/07_tax/sum_${i}.txt
    done
    sed -i 's/(//g;s/)//g;s/\"//g;s/\#//g;s/\/Chloroplast//g' result/07_tax/sum_*.txt


## 1.12. Space clearance and data submission
    rm -rf temp/*.fq
    parallel -j 8 --xapply \
      "gzip {1}" \
      ::: seq/*

    # Calculate the MD5 hash on both ends for data submission
    cd seq
    md5sum *.fq.gz > ../result/md5sum.txt
    cat ../result/md5sum.txt
    md5sum -c ../result/md5sum.txt
    cd ..



# 2. Diversity and Species Analysis in R


## 2.1. Alpha Diversity

### 2.1.1 Alpha diversity box plot
    for j in `head -n 1 result/metadata.tsv |cut -f 2,3,7`;do
      for i in `head -n1 result/05_alpha/vegan.txt|cut -f 2-`;do
        Rscript ~/scripts/alpha_boxplot.R --alpha_index ${i} \
          --input result/05_alpha/vegan.txt --design result/metadata.tsv \
          --group ${j} --output result/05_alpha/ \
          --width 89 --height 59
      done
    done
    mv alpha_boxplot_TukeyHSD.txt result/05_alpha/

### 2.1.2 Rarefaction Curve
    Rscript ~/scripts/alpha_rare_curve.R \
      --input result/05_alpha/alpha_rare.txt --design result/metadata.tsv \
      --group Fermenter --output result/05_alpha/ \
      --width 178 --height 118

### 2.1.3 Diversity Venn diagram
    bash ~/scripts/sp_vennDiagram.sh \
      -f result/05_alpha/otu_group_exist2.txt \
      -a FFS1 -b FFS2 -c FFS3 \
      -w 3 -u 3 \
       -p FFS1_FFS2_FFS3
    
## 2.2. Beta Diversity
    for j in `head -n 1 result/metadata.tsv |cut -f 2,3,7`;do
      for i in bray_curtis euclidean unifrac jaccard;do
        Rscript ~/scripts/beta_cpcoa.R \
          --input result/06_beta/${i}.txt --design result/metadata.tsv \
          --group ${j} --output result/06_beta/${i}.${j}.cpcoa.pdf \
          --width 89 --height 59
      done
    done
      
    
## 2.3. Taxonomy

### 2.3.1 Stackplot
    Rscript ~/scripts/tax_stackplot.R \
      --input result/07_tax/sum_g.txt --design result/metadata.tsv \
      --group Stage --output result/07_tax/sum_g.stackplot \
      --legend 10 --width 178 --height 118

    for i in p c o f g s; do
      Rscript ~/scripts/tax_stackplot.R \
        --input result/07_tax/sum_${i}.txt --design result/metadata.tsv \
        --group Stage --output result/07_tax/sum_${i}.stackplot \
        --legend 10 --width 118 --height 88
    done
        
    for i in p c o f g s; do
      Rscript ~/scripts/mean.R --input result/07_tax/sum_${i}.txt \
        --design result/metadata.tsv \
        --group Stage --thre 0 \
        --output result/07_tax/sum_${i}_mean.txt
      tail -n+2 result/07_tax/sum_${i}_mean.txt |sort -k2,2nr |sed "1i `sed -n 1p result/07_tax/sum_${i}_mean.txt`" \
        > result/07_tax/sum_${i}_mean_sort.txt
      rm result/07_tax/sum_${i}_mean.txt
    done
      
    

### 2.3.2 circlize
    for i in p c o f g s; do
    Rscript ~/scripts/tax_circlize.R \
      --input result/07_tax/sum_${i}.txt --design result/metadata.tsv \
      --group Stage --legend 10
    mv circlize.pdf result/07_tax/sum_${i}.circlize.pdf
    mv circlize_legend.pdf result/07_tax/sum_${i}.circlize_legend.pdf
    done


# 3. Comparison of Differences

## 3.1. Differential Analysis in R
### 3.1.1 Difference comparison
    mkdir -p result/09_compare/
    compare="BFS1-BFS2"
    Rscript ~/scripts/compare.R \
      --input result/otutab.txt --design result/metadata.tsv \
      --group Stage --compare ${compare} --threshold 0.1 \
      --method wilcox --pvalue 0.05 --fdr 0.05 \
      --output result/09_compare/

### 3.1.2 Volcano diagram
    Rscript ~/scripts/compare_volcano.R \
      --input result/09_compare/${compare}.txt \
      --output result/09_compare/${compare}.volcano.pdf \
      --width 89 --height 59

### 3.1.3 Map of Manhattan
    bash ~/scripts/compare_manhattan.sh -i result/09_compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/07_tax/sum_p.txt \
       -w 183 -v 59 -s 7 -l 10 \
       -o result/09_compare/${compare}_p

   bash ~/scripts/compare_manhattan.sh -i result/09_compare/${compare}.txt \
       -t result/taxonomy.txt \
       -p result/07_tax/sum_g.txt \
       -w 450 -v 150 -s 10 -l 30 -L Genus \
       -o result/09_compare/${compare}_g
       
    rm Rplots.pdf


## 3.2. LEfSe
    Rscript ~/scripts/format2lefse.R --input result/otutab.txt \
      --taxonomy result/taxonomy.txt --design result/metadata.tsv \
      --group Stage --threshold 0.1 \
      --output result/09_compare/LEfSe
    
    mkdir -p ${wd}/result/09_compare/LEfSe && cd ${wd}/result/09_compare/LEfSe
    cp ${wd}/result/09_compare/LEfSe.txt ./
    
    conda activate meta_lefse
    lefse-format_input.py LEfSe.txt input.in -c 1 -o 1000000
    run_lefse.py input.in input.res -l 2
    cat input.res |grep '^[A-Z][a-z]*.[A-Z][a-z]*.[A-Z][a-z]*.[A-Z][a-z]*.[A-Z][a-z]*.[A-Z][a-z]*' > input_genus.res
    # lefse-plot_cladogram.py input_genus.res cladogram_genus.pdf --format pdf
    
    # Create a bar chart showing all differences
    lefse-plot_res.py input.res res.pdf  --subclades -5 --format pdf
    lefse-plot_res.py input_genus.res res_genus.pdf  --subclades -5 --format pdf
    
    # Plot a bar chart for individual features
    for i in `cat input_genus.res |cut -f1`;do
      lefse-plot_features.py -f one --feature_name ${i} \
         --format pdf input.in input.res ${i}.pdf
    done
