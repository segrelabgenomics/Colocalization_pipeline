base_dir/eCaviar_table/

*Step1:* Generate Table 1 with QC

* Parse through eCAVIAR results
    * run_ecaviar_results_table.sh
        Required inputs:
        ``--ecaviar_directory``: Directory consisting of eCAVIAR outputs
        ``--gwas_trait``: GWAS trait of interest
        ``--output_file``: Name of output file
        ``--incorrect_num_file``:Name of output file consisting of loci with incorrect matching of col and z-score files
        
    * run_curate_ecaviar_table_pretable1.sh 
        Required inputs:
        ``--ecaviar_table``: eCAVIAR table with all aggregate results
        ``--tissue``: Tissue
        ``--qtl_type``: QTL Type
        ``--gene``: Gene 
        ``--gwas_trait``: GWAS trait of interest
        ``--output_directory``: Name of output file
        
        - Required file: /base_dir/eCaviar_table/tables/MANIFEST_gwas_trait_slope_se.txt
        This file contains tissue;QTL;gene_id (gene or splice event) without header
        
    * run_concatenate_output_tables.sh
        -Create a directory:
        /base_dir/eCaviar_table/tables/Intermediate_files/results_concatenated 
        /base_dir/eCaviar_table/tables/Intermediate_files/eqtl 
        /base_dir/eCaviar_table/tables/Intermediate_files/sqtl 
        Required inputs:
        ``--base_dir_eqtl``: Directory containing outputs for eqtls
        ``--base_dir_sqtl``: Directory containing outputs for sqtls
        ``--output_directory``: Directory to place parsed files
        ``--output_filename``: Name of output file
        
        -Once the script has successfully run. Move files to sQTL and eQTL folder respectively
    
   * run_curate_table1.sh
        Required inputs:
        ``--ecaviar_output_table``: eCaviar output file which is concatenated
        ``--gwas_file``: GWAS file
        ``--output_filename``: Name of output file
        
    * Generate Table 1
        finish_table_1_parsing.py: This requires multiple input files set up within the script. User needs to change everytime the script is exceuted.
        
    * run_add_significat_column.sh 
        Required inputs:
        ``--cutoff_file``: eCaviar Table1
        ``--outputs_dir``: Name of output directory
        
    * change_qc_cutoff_table1.py - adapt to your desired cut-off
    * change_final_header_table1.py - standardize headers
    

*Step 2*: Generate a summary table

   * Generate Table 2
        finish_table_2_parsing.py: This requires multiple input files set up within the script. User needs to change everytime the script is exceuted.
