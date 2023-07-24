from os import path, makedirs, rename, remove

# Show the rulegraph
#  snakemake -s Snakefile --rulegraph | dot -Tpng > rulegraph.png

configfile: "config/config.yaml"

onsuccess:
    print("Workflow completed without any error")

onerror:
    print("An error occurred")

# Should an excel file with barcodes and sample names be coverted to TSV ?
# Please note that it is assumed that the 5th and 2nd columns of sheet one
#  in the Excel file are the sample names and barcodes, respectively
isExcel=config['isExcel']  # True or False
EXCEL_FILE=config['EXCEL_FILE']

SAMPLES=config['SAMPLES']

TAXON_LEVELS={'p':'Phylum', 'c': 'Class', 'o': 'Order','f': 'Family','g': 'Genus', 's': 'Species'}

rule all:
    input: 
        "08.Quality_check/qual.txt",
        "16.Alpha_diversity/alpha.txt",
        "18.Beta_diversity/beta_diversity.done",
        expand("20.Summarize_taxonomy/{taxon_level}_table.txt", taxon_level=TAXON_LEVELS.values())
      

rule Join_fastq:
    input:
        forward_index = config["INDEX1"], #"01.raw_data/index1.fastq",
        reverse_index = config["INDEX2"] #"01.raw_data/index2.fastq"
    output: "02.Join_fastq/indices.fastq"
    message: "Merging the forward ({input.forward_index}) and reverse ({input.reverse_index}) indices"
    params:
        program=config['programs_path']['usearch']
    threads: 1
    shell:
        """
        {params.program} \
                 -fastq_join {input.forward_index} \
                 -reverse {input.reverse_index} \
                 -join_padgap "" \
                 -threads {threads} \
                 -fastqout {output}
        """

# Convert barcodes file from TSV to FASTA
rule Parse_barcodes:
    input: EXCEL_FILE #if isExcel else "01.raw_data/sample2barcode.tsv"
    output: "03.Parse_barcodes/bar.fasta"
    message: "Converting the barcodes file from TSV to FASTA.." 
    params:
        isExcel=isExcel,
        program=config['programs_path']['xlsx2csv']
    threads: 1
    shell:
        """
        ISEXCEL={params.isExcel}
        if [ ${{ISEXCEL}} == True ]; then

             {params.program} -d "tab" -s 1 {input} | \
              awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print($5,$2)}}' | \
              awk '{{print ">"$1"\\n"$2}}' > {output}

        else

            awk '{{print ">"$1"\\n"$2}}' {input} > {output}

        fi
        """
        

rule Reformat_barcodes:
    input: rules.Parse_barcodes.output
    output: "04.Reformat_barcodes/bar.fasta"
    message: "Formatting the barcode sequences to the orientation required by usearch...."
    threads: 1
    script:
        "./reformat_barcode.py"

# Demultiplex the reads per sample using the reformated barcodes
rule Demultiplex:
    input: 
        barcodes=rules.Reformat_barcodes.output,
        indices=rules.Join_fastq.output,
        forward=config['READ1'], # "01.raw_data/read1.fastq",
        rev=config['READ2'] #"01.raw_data/read2.fastq"
    output:
        forward="05.Demultiplex/demux_R1.fastq",
        rev="05.Demultiplex/demux_R2.fastq"
    message: "Demultiplexing.... Retreiving analysis specific samples from the entire sequence run..."
    params:
        program=config['programs_path']['usearch'],
        out_dir=lambda w, output: path.dirname(output[0])
    threads: 1
    shell:
        """       
        # Demultiplex
        {params.program} \
                -fastx_demux {input.forward} \
                -reverse {input.rev} \
                -index {input.indices} \
                -barcodes {input.barcodes} \
                -fastqout {output.forward} \
                -output2 {output.rev}
        """

rule Trim_primers:
    input:
        forward=rules.Demultiplex.output.forward,
        rev=rules.Demultiplex.output.rev
    output:
        forward="06.Trim_primers/reads_R1.fastq",
        rev="06.Trim_primers/reads_R2.fastq"
    message: "Trimming out primer sequences..."
    params:
        program=config['programs_path']['usearch'],
        len2trim=config['parameters']['fastx_truncate']['len2trim'] # primer length
    threads: 1
    shell:
        """
        # Trim foward and reverse primers
        {params.program} -fastx_truncate \
          {input.forward} \
          -stripleft {params.len2trim} \
          -fastqout {output.forward} &&
        {params.program} -fastx_truncate \
           {input.rev} \
          -stripleft {params.len2trim} \
          -fastqout {output.rev}
        """

rule Merge_pairs:
    input: 
        forward=rules.Trim_primers.output.forward,
        rev=rules.Trim_primers.output.rev
    output: "07.Merge_pairs/merged.fastq"
    message: "Merging pairs..."
    log: "logs/Merge_pairs/merged.txt"
    params:
        program=config['programs_path']['usearch'],
        max_diff=config['parameters']['fastq_mergepairs']['max_diff'],
        trunc_tail=config['parameters']['fastq_mergepairs']['trunc_tail'],
        min_overlap_len=config['parameters']['fastq_mergepairs']['min_overlap_len'],
        min_merge_len=config['parameters']['fastq_mergepairs']['min_merge_len']
    shell:
        """
        # Merge pairs
        {params.program} -fastq_mergepairs \
           {input.forward} -reverse {input.rev} \
           -fastq_maxdiffs {params.max_diff} \
           -fastq_trunctail {params.trunc_tail} \
           -fastq_minovlen {params.min_overlap_len} \
           -fastq_minmergelen {params.min_merge_len} \
           -fastqout {output} > {log} 2>&1
        """

rule Quality_check:
    input: rules.Merge_pairs.output
    output: "08.Quality_check/qual.txt"
    message: "Just quality checking... See the output here {output}"
    log: "logs/Quality_check/Quality_check.txt"
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        {params.program} -fastq_eestats2 \
        {input} -output {output} > {log} 2>&1
        """

rule Filter_seqs:
    input: rules.Merge_pairs.output
    output: "09.Filter_seqs/Filtered.fa"
    message: "Filtering and trimming away low quality sequences..."
    log: "logs/Filter_seqs/Filter_seqs.txt"
    params:
        program=config['programs_path']['usearch'],
        max_ee=config['parameters']['fastq_filter']['max_ee'],
        trunc_len=config['parameters']['fastq_filter']['trunc_len']
    shell:
        """
        # Quality filter & trim to length 350 (chosen after review of qual report)
        {params.program} -fastq_filter {input} \
            -fastq_maxee {params.max_ee} \
            -fastq_trunclen {params.trunc_len} \
            -fastq_maxns 0 \
            -relabel Filt \
            -fastaout {output} > {log} 2>&1
        """

rule Dereplicate_seqs:
    input: rules.Filter_seqs.output
    output: "10.Dereplicate_seqs/uniques.fa"
    message: "Dereplicating sequences...."
    log: "logs/Dereplicate_seqs/Dereplicate_seqs.txt"
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        # Find unique read sequences and their abundance
        {params.program} -fastx_uniques {input} \
           -sizeout \
           -relabel Uniq \
           -fastaout {output} > {log} 2>&1
        """

rule Remove_singletons:
    input: rules.Dereplicate_seqs.output
    output: "11.Remove_singletons/uniques.sorted.min2.fa"
    message: "Removing singletons..."
    log: "logs/Remove_singletons/Remove_singletons.txt"
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        # Removing singletons i.e. sequences that appear only once
        {params.program} -sortbysize {input} \
           -fastaout {output} -minsize 2 > {log} 2>&1
        """

rule Denoise:
    input: rules.Remove_singletons.output
    output: 
        features="12.Denoise/features.fa",
        table="12.Denoise/unoise3.txt"
    message: "Denoising...."
    log: "logs/Denoise/Denoise.txt"
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        {params.program} -unoise3 {input} \
            -zotus {output.features} \
            -tabbedout {output.table} > {log} 2>&1
        """

rule Make_ASV_table:
    input: 
        features=rules.Denoise.output.features,
        fastq=rules.Merge_pairs.output
    output: "13.Make_ASV_table/feature_table.txt"
    message: "Making the ASV table..."
    log: "logs/Make_ASV_table/Make_ASV_table.txt"
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        {params.program} -otutab {input.fastq} \
           -otus {input.features} \
           -otutabout {output}  > {log} 2>&1
        """

# Get the refaction depth ike so:
# awk 'BEGIN{FS=OFS="\t"} NR>1{a+=$2; b+=$3; c+=$4; d+=$5; e+=$6; f+=$7; g+=$8; h+=$9; i+=$10; j+=$11} END{printf "%1.f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\t%1.f\n", a,b,c,d,e,f,g,h,i,j}' 13.Make_ASV_table/feature_table.txt

rule Summarize_table:
    input: rules.Make_ASV_table.output
    output: "14.Summarize_table/summary.txt"
    message: "Generating a summary file {output} for the ASV table. Check this file out to get the number for rarefying to even depth."
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        {params.program} -otutab_stats {input} -output {output}
        """

rule Rarefy:
    input: 
        table=rules.Make_ASV_table.output,
        summary=rules.Summarize_table.output
    output: "15.Rarefy/feature_table.txt"
    message: "Rarefying the ASV table to even depth at {params.depth}..."
    params:
        program=config['programs_path']['usearch'],
        depth=config['parameters']['otutab_rare']['depth']
    shell:
        """
        {params.program}  -otutab_rare {input.table} \
           -sample_size {params.depth} \
           -output {output}
        """

rule Alpha_diversity:
    input: rules.Rarefy.output
    output: "16.Alpha_diversity/alpha.txt"
    message: "Performing Alpha diversity analysis..."
    params:
        program=config['programs_path']['usearch']
    shell:
        "{params.program} -alpha_div {input} -output {output}" 

rule Make_tree:
    input: rules.Denoise.output.features
    output: "17.Make_tree/features.tree"
    message: "Generating an ASV tree..."
    params:
        program=config['programs_path']['usearch'],
        outdir=lambda w,output: path.dirname(output[0])
    shell:
        """
        {params.program} -calc_distmx {input} \
             -tabbedout {params.outdir}/distance_matrix.txt && \
        {params.program} -cluster_aggd  {params.outdir}/distance_matrix.txt \
             -treeout {output}
        """

rule Beta_diversity:
    input: table=rules.Rarefy.output, tree=rules.Make_tree.output
    output:  
        directory("18.Beta_diversity/"), 
        touch("18.Beta_diversity/beta_diversity.done") 
    message: "Performing Beta diversity analysis..."
    params:
        program=config['programs_path']['usearch']
    shell:
        """
        {params.program} -beta_div {input.table} \
            -tree {input.tree} \
            -filename_prefix {output[0]}/ \
            -metrics jaccard,bray_curtis,unifrac
        """

rule Assign_taxonomy:
    input: rules.Denoise.output.features
    output: "19.Assign_taxonomy/sintax.txt"
    message: "Assigning taxonomy to the representative ASV sequences..."
    log: "logs/Assign_taxonomy/sintax.log"
    params:
        program=config['programs_path']['usearch'],
        database=config['parameters']['sintax']['database'],
        cutoff=config['parameters']['sintax']['cutoff']
    shell:
        """
        {params.program} -sintax {input} \
            -db {params.database} -strand both  \
            -sintax_cutoff {params.cutoff} -tabbedout {output}
        """

#TAXON_LEVELS={'p':'Phylum', 'c': 'Class', 'o': 'Order','f': 'Family','g': 'Genus', 's': 'Species'}
rule Summarize_taxonomy:
    input: 
        table=rules.Make_ASV_table.output,
        taxonomy=rules.Assign_taxonomy.output
    output: expand("20.Summarize_taxonomy/{taxon_level}_table.txt", taxon_level=TAXON_LEVELS.values())
    message: "Performing Beta diversity analysis..."
    params:
        program=config['programs_path']['usearch']
    run:
        for rank, name in TAXON_LEVELS.items():
            shell("""
                  {params.program} -sintax_summary {input.taxonomy} \
                      -otutabin {input.table} \
                      -rank {rank} \
                      -output 20.Summarize_taxonomy/{name}_table.txt
                  """)
    

