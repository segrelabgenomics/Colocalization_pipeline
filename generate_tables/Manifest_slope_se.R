args <- commandArgs(trailingOnly = TRUE)
trait <- ifelse(length(args) >= 1, args[1], Sys.getenv("POAG_EUR", "TRAIT_NAME"))
input_file <- ifelse(length(args) >= 2, args[2], Sys.getenv("ECAVIAR_RESULTS_TABLE", file.path("generate_tables/tables", paste0(trait, "_eCAVIAR_results_table.txt"))))
output_file <- ifelse(length(args) >= 3, args[3], Sys.getenv("SLOPE_SE_MANIFEST", file.path("generate_tables/tables", paste0("MANIFEST_", trait, "_slope_se.txt"))))

manifest_slope_se <- read.delim2(input_file)
manifest_slope_se <- manifest_slope_se[,c("tissue","qtl","gene_id")]
manifest_slope_se <- paste0(manifest_slope_se$tissue,";",manifest_slope_se$qtl,";",manifest_slope_se$gene_id)
write.table(manifest_slope_se,
            output_file,
            row.names = F,
            col.names = F,
            quote = F,
            append = F)
