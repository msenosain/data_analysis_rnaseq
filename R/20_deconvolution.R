# Deconvolution

# MCP counter
mcp_dcv <- MCPcounter::MCPcounter.estimate(tpm_all, featuresType = "HUGO_symbols")
write.table(mcp_dcv, file = "mcp_dcv.txt", sep = "\t",
            row.names = TRUE, quote = FALSE)


# quanTIseq
qts_dcv <- immunedeconv::deconvolute(tpm_all, "quantiseq")
write.table(qts_dcv, file = "qts_dcv.txt", sep = "\t",
            row.names = FALSE, quote = FALSE)


# CIBERSORT and xCell done online