.PHONY: run clean

run:
	snakemake -pr --keep-going --rerun-incomplete --cores 10 

clean:
	rm -rf 02.Join_fastq/  03.Parse_barcodes/ 04.Reformat_barcodes/ 05.Demultiplex/ 06.Trim_primers/ 07.Merge_pairs/ 08.Quality_check/ 09.Filter_seqs/ 10.Dereplicate_seqs/ 11.Remove_singletons/ 12.Denoise/ 13.Make_ASV_table/ 14.Summarize_table/ 15.Rarefy/ 16.Alpha_diversity/ 17.Make_tree/ 18.Beta_diversity/ 19.Assign_taxonomy/ 20.Summarize_taxonomy/  
