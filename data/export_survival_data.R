# Export the public datasets used in the applications from the R package survival.
# Run from the data/ folder:  Rscript export_survival_data.R
# Then run:                   python prepare_data.py
# The analyses were run with survival 3.8-3 (R 4.4.3).
suppressPackageStartupMessages(library(survival))

rotterdam <- data.frame(rownames = rownames(survival::rotterdam), survival::rotterdam, check.names = FALSE)
gbsg <- data.frame(rownames = rownames(survival::gbsg), survival::gbsg, check.names = FALSE)
write.csv(rotterdam, "rotterdam.csv", row.names = FALSE, quote = FALSE, na = "")
write.csv(gbsg, "gbsg.csv", row.names = FALSE, quote = FALSE, na = "")
write.csv(survival::nwtco, "nwtco.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "R_SESSION_INFO.txt")
writeLines(capture.output(citation("survival")), "SURVIVAL_CITATION.txt")
cat("survival", as.character(packageVersion("survival")), "\n")
cat("rotterdam", nrow(rotterdam), "gbsg", nrow(gbsg), "nwtco", nrow(survival::nwtco), "\n")
