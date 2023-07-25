.PHONY: usage run rarefy_depth clean

usage:
	@echo "make rarefy_depth # Get the rarefaction depth"
	@echo "make run # Run the complete pipeline"
	@echo "make clean # Delete all the generated directories"

run:
	@echo "Run the complete pipeline"
	snakemake -pr --keep-going --rerun-incomplete --cores 10 

rarefy_depth:
	@echo "Get the refaction depth for divesity analysis by viewing the file '14.Summarize_table/summary.txt' "
	snakemake -pr --keep-going --rerun-incomplete --cores 10  "14.Summarize_table/summary.txt"

clean:
	@echo "Removing all the generated ouptuts"
	rm -rf 02.Join_fastq/  03.Parse_barcodes/ 04.Reformat_barcodes/ 05.Demultiplex/ 06.Trim_primers/ 07.Merge_pairs/ 
	rm -rf 08.Quality_check/ 09.Filter_seqs/ 10.Dereplicate_seqs/ 11.Remove_singletons/ 12.Denoise/ 13.Make_ASV_table/ 
	rm -rf 14.Summarize_table/ 15.Rarefy/ 16.Alpha_diversity/ 17.Make_tree/ 18.Beta_diversity/ 19.Assign_taxonomy/ 20.Summarize_taxonomy/  
