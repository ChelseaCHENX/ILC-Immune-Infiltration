# https://github.com/omnideconv/immunedeconv
# consensusTME, https://aacrjournals.org/cancerres/article/79/24/6238/639705/Comprehensive-Benchmarking-and-Integration-of
# Finally, certain methods can be used with custom signatures, consisting of either a signature matrix or signature genes for the cell types of interest. 

# source activate deconvolution; R
rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'
setwd(rootdir)
source('codes/functions_brief.R')

#---------------------------------- purity
purity = immunedeconv::deconvolute_estimate(gex)
write.csv(purity,'data/deconv/deconvolute_estimate.csv', quote=F)

#---------------------------------- deconvolution scores
# gex = read.csv('data/gex/tpm.csv', row.names=1)
ma_sigs = parse_gmt('data/geneset/Macrophages.gmt')

# metabric
fpath = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/yap/data/normcts/metabric.nolan.expr.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_metabric.csv', quote=F)

# scan-b
fpath = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/yap/data/normcts/scanb_log2cpm.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_scanb.csv', quote=F)

# poetic
fpath = 'data/clinical/peotic_gex.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_poetic.csv', quote=F)

sessionInfo()
