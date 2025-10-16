library(CellChat)
library(tictoc)
library(future)
library(tidyverse)
# future::plan("multisession", workers = 20)
# options(future.globals.maxSize = 30 * 1024 ^ 3) # for 50 Gb RAM
# future::plan("sequential")
dir_res <- '/volumes/USR1/yyan/project/tnbc_pre_atlas/rds_rna-integrate/pat102/ecotype/cellchat'

cellchat <- read_rds(file.path(dir_res, 'cellchat.rds'))
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.human

cellstatenames <- levels(cellchat@idents)
#------------------- ~~~ Preprocessing ~~~ -------------------  
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
write_rds(cellchat, file.path(dir_res, 'cellchat.rds'))

cellchat <- identifyOverExpressedInteractions(cellchat)
write_rds(cellchat, file.path(dir_res, 'cellchat.rds'))
if (F) {
  ## // Did not intend to use the PPI-projected values
  cellchat <- projectData(cellchat, PPI.human)
}
#------------------- ~~~ Infering cell-cell interaction ~~~ -------------------  

# // Do use the PPI-projected values 
# // Adjust for the abundance of cell states
tic('computeCommunProb...')
cellchat <- computeCommunProb(cellchat, 
                              raw.use = T, 
                              population.size = TRUE)
toc()
write_rds(cellchat, file.path(dir_res, 'cellchat.rds'))
write_lines('computeCommunProb', file.path(dir_res, 'log.txt'))

write_rds(cellchat, file.path(dir_res, 'cellchat.0.rds')) # log


cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
write_rds(cellchat, file.path(dir_res, 'cellchat.rds'))
write_lines('computeCommunProbPathway', file.path(dir_res, 'log.txt'), append = T) # log


cellchat <- read_rds(file.path(dir_res, 'cellchat.rds'))