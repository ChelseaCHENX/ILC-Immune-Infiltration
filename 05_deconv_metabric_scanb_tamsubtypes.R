rootdir = '/Users/chelseachen/Library/CloudStorage/OneDrive-UniversityofPittsburgh/01_projects/ILCimmune'
setwd(rootdir)
source('codes/functions_brief.R')

gex = read.csv('data/gex/tpm.csv', row.names=1)
#---------------------------------- purity
purity = immunedeconv::deconvolute_estimate(gex)
write.csv(purity,'data/deconv/deconvolute_estimate.csv', quote=F)

#---------------------------------- deconvolution scores
ma_sigs = parse_gmt('data/geneset/Macrophages.gmt')

# metabric
fpath = 'data/public/metabric.nolan.expr.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_metabric.csv', quote=F)

# scan-b
fpath = 'data/public/scanb_log2cpm.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_scanb.csv', quote=F)

# poetic
fpath = 'data/public/peotic_gex.csv'
gex <- fread(fpath, data.table = FALSE)
row.names(gex) <- gex[, 1]
gex <- gex[, -1, drop = FALSE]

res = deconvolute_consensus_tme_custom(as.matrix(gex), ma_sigs, 
stat_method='ssgsea')
write.csv(res,'data/deconv/deconvolute_tams_poetic.csv', quote=F)

# sessionInfo()

