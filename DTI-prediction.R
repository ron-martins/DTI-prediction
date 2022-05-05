######################################################################
######################################################################
rm(list = ls())

library(ggplot2)
library(dplyr)
library(corrplot)
library(RColorBrewer)
library(seqinr)
library(Biostrings)
library(data.table)
library(igraph)
library(pracma)
library(glue)
library(ggpubr)
library(ggsci)
library(splitstackshape)

######################################################################
######################################################################


######################################################################
############ Database ############

### BindingDB_All database -> importing ###
## ## Due to the size of the file, it is necessary to download it directly from the BindingDB database.
ref.all <- read.csv('BindingDB_All.tsv', '\t', stringsAsFactors = F, header = T)

### Selecting columns and saving the Referenced database ###
refdb <- ref.all %>%
  select(c("BindingDB.Reactant_set_id","Ligand.SMILES","Ligand.InChI.Key","BindingDB.Ligand.Name",
           "Target.Name.Assigned.by.Curator.or.DataSource","Target.Source.Organism.According.to.Curator.or.DataSource",
           "Ki..nM.","Kd..nM.","IC50..nM.","EC50..nM.",
           "Curation.DataSource","PDB.ID.s..for.Ligand.Target.Complex",
           "PubChem.CID","PubChem.SID","ChEBI.ID.of.Ligand","ChEMBL.ID.of.Ligand",
           "DrugBank.ID.of.Ligand","IUPHAR_GRAC.ID.of.Ligand","KEGG.ID.of.Ligand","ZINC.ID.of.Ligand",
           "UniProt..SwissProt..Recommended.Name.of.Target.Chain","UniProt..SwissProt..Primary.ID.of.Target.Chain", "UniProt..SwissProt..Entry.Name.of.Target.Chain"))


### Tidying, applying filters and fixing issues ### 
refdb <- refdb %>%
  filter(Target.Source.Organism.According.to.Curator.or.DataSource =='Homo sapiens') %>%
  transform(Ki..nM. = as.numeric(as.character(Ki..nM.)), Kd..nM. = as.numeric(as.character(Kd..nM.)),
            IC50..nM. = as.numeric(as.character(IC50..nM.)), EC50..nM. = as.numeric(as.character(EC50..nM.))) %>%
  filter(Ki..nM. < 10000 | Kd..nM. <10000 | IC50..nM. < 10000 | EC50..nM. < 10000) %>%
  rename(InChIKey = Ligand.InChI.Key, SMILES = Ligand.SMILES, UniProt = UniProt..SwissProt..Primary.ID.of.Target.Chain) %>%
  filter(!is.na(UniProt) & !is.na(InChIKey) & !is.na(SMILES)) %>%
  filter(!(UniProt == "") & !(InChIKey == "") & !(SMILES == ""))


######################################################################
############ natural product databases (NuBBE) ############



### NuBBE list (plant) -> importing ###
nubbedb <- read.csv('../data/plants_NuBBE.inchikey', header = F, stringsAsFactors = F)
names(nubbedb) <- 'key' 

### Comparing NuBBE list to Referenced databases (filtered) ###
### Selecting the intersting columns ###
nubbe.ref <- refdb[refdb$InChIKey %in% nubbedb$key, ] %>%
  select(c(BindingDB.Reactant_set_id, BindingDB.Ligand.Name, InChIKey, SMILES, 
           DrugBank.ID.of.Ligand, PubChem.CID,
           UniProt, UniProt..SwissProt..Entry.Name.of.Target.Chain, UniProt..SwissProt..Recommended.Name.of.Target.Chain))


######################################################################
############ selecting targets from UniProt database  ############

### importing and tidying UniProt database ###
UniProt.df <- data.frame(names(readAAStringSet('uniprot_sprot.fasta')), 
                         paste(readAAStringSet('uniprot_sprot.fasta')))

colnames(UniProt.df) <- c('seq_name', 'sequence')

UniProt.seq <- unique(subset(UniProt.df %>%
                               cSplit(splitCols = "seq_name",
                                      sep = "|", direction = "long") %>%
                               filter(!is.na(seq_name)), seq_name!= "sp") %>% 
                        distinct(sequence, .keep_all = T))

UniProt.seq[] <- lapply(UniProt.seq, as.character)

### selecting UniProt code from NuBBE list ###
target.seq = matrix(data = NA, nrow = nrow(target.list), ncol = 1, dimnames = c(target.list, "Sequence"))
target.list = unique(nubbe.ref %>% select(c(UniProt)))

for(i in 1:nrow(target.seq)){
  x = target.list[i, ]
  target.seq[i, ] = filter(UniProt.seq, seq_name == x)$sequence
}


######################################################################
############ importing and/or selecting other datasets ############

### Drug list (Inchikey) ###
nubbe.ref.drug <- unique(nubbe.ref %>% select(InChIKey))
names(nubbe.ref.drug) <- 'key'

### Drug Fingerprint from PaDEL- PubChem software ###
FP.pub <- read.csv('../data/NuBBE_ref_drugs-Pub.csv')
FP.pub <- cbind(nubbe.ref.drug, FP.pub) %>% select(-Name)


######################################################################
############ all functions ############ 

normav <- function(x) 
{
  return ((x - min(x)) / (max(x) - min(x)))
}

Jaccard = function (x, y) 
{
  M.11 = sum(x == 1 & y == 1)
  M.10 = sum(x == 1 & y == 0)
  M.01 = sum(x == 0 & y == 1)
  return (M.11 / (M.10 + M.01 + M.11))
}

NBI <- function(list, FP.pub, target.seq, ll, ld, lt)
{
  
  ## list -> dataset input
  ## FP.pub -> fingerprint from PaDEL
  ## target.seq -> target sequence FASTA (from UniProt)
  ## ll, ld, lt -> CMF parameters
  
  g <- graph_from_data_frame(unique(list %>% select(c("InChIKey", "UniProt"))),
                             directed = F)
  V(g)$type <- bipartite.mapping(g)$type
  
  V(g)$K <- ego_size(g, order = 1, nodes = V(g), mindist = 1)
  V(g)$K1 <- 1/(ego_size(g, order = 1, nodes = V(g), mindist = 1))
  g.dsub <- induced.subgraph(graph = g, vids = V(g)$type == F) 
  g.tsub <- induced.subgraph(graph = g, vids = V(g)$type == T) 
  
  W <- matrix(ncol = vcount(g.dsub), nrow = vcount(g.dsub)) 
  for(p in 1:dim(W)[1]){
    for(q in 1:dim(W)[2]){
      W[p,q] =  V(g)[q]$K1*(sum(intersection(neighbors(g, p), neighbors(g, q))$K1))
    }  
  }
  
  WF <- matrix(ncol = vcount(g.tsub), nrow = vcount(g.dsub), 
               dimnames = list(V(g.dsub)[]$name, WFcnames <- V(g.tsub)[]$name)) 
  for(d in 1:nrow(WF)){
    for(t in 1:ncol(WF)){
      WF[d,t] =  (W %*% g[V(g)$type == F, V(g)$type == T][,t])[d]
    }  
  }
  
  for(p in 1:ncol(WF)){
    WF[ ,p] = normav(WF[ ,p]) 
  }
  WF <- reshape2::melt(as.matrix(WF[hclust(dist(WF))$order, hclust(dist(t(WF)))$order])) %>%
    setNames(c("Drug_ID", "Target_ID", "value")) %>%
    cbind(method = 'NBI')
  
  return(WF)
}

Simil.d <- function (list, FP)
{
  ## list -> dataset input
  ## FP -> fingerprint from PaDEL
  
  #transpose:
  FP.t <- data.frame(t(FP))
  names(FP.t) <- as.matrix(FP.t[1, ])
  FP.t <- FP.t[-1, ]
  
  # cleaning the transpose with ONLY drugs from dataset:
  FP.t <- FP.t[, (colnames(FP.t) %in% unique(list$InChIKey))]
  
  #output matrix:
  m = matrix(data = NA, nrow = length(FP.t), ncol = length(FP.t))
  
  #applying Jaccard function: 
  for (r in 1:length(FP.t)) {
    for (c in 1:length(FP.t)) {
      if (c == r) {
        m[r,c] = 1
      } else {
        m[r,c] = Jaccard(FP.t[,r], FP.t[,c])
      }
    }
  }
  
  #naming the similarity drugxdrug matrix: 
  for(i in as.data.frame(colnames(FP.t))){
    colnames(m)=i
    rownames(m)=i
  }
  
  return(m)
  
}

Simil.t <- function (list, target.seq)
{
  ## list -> dataset input
  ## target.seq -> target sequence FASTA (from UniProt)
  
  #sequence matrix:
  target.seq = matrix(data = NA, nrow = nrow(list), ncol = 1, dimnames = c(list, "Sequence"))
  
  #obtain sequences from SwissProt databank:
  choosebank("swissprot")
  for(i in 1:nrow(target.seq)){
    x = paste("AC=", list[i, ], sep='')
    target.seq[i, ] = c2s(getSequence(query(x)$req[[1]]))  }
  
  # cleaning the sequences with ONLY targets from dataset:
  #target.seq$UniProt <- row.names(target.seq)
  #target.seq <- target.seq[target.seq$UniProt %in% unique(list$UniProt), ]
  
  #aux matrix:
  target.simil = matrix(data = NA, nrow = nrow(target.seq), 
                        ncol = nrow(target.seq), 
                        dimnames = c(target.seq %>% select(-Sequence), 
                                     target.seq %>% select(-Sequence)))
  
  #applying similarity:
  data("BLOSUM62")
  for(a in 1:nrow(target.simil)){
    for(b in 1:ncol(target.simil))
    {
      target.simil[a,b] = pairwiseAlignment(target.seq[a, 1], 
                                            target.seq[b, 1], 
                                            type = 'local', 
                                            substitutionMatrix=BLOSUM62, scoreOnly = T)
    }  
  }
  
  #output matrix:
  target.simil.norm <- target.simil
  
  #normalizing the similarity scores (0 to 1)
  for(p in 1:ncol(target.simil)){
    target.simil.norm[ ,p] = normav(target.simil[ ,p]) 
  }
  
  return(target.simil.norm)
  
}

CMF <- function(list, FP.pub, target.seq, ll, ld, lt)
{
  ## list -> dataset input
  ## FP.pub -> fingerprint from PaDEL
  ## target.seq -> target sequence FASTA (from UniProt)
  ## ll, ld, lt -> CMF parameters 
  
  #sleecting columns from dataset:
  Y <- as.matrix(table(unique(list %>% select(c("InChIKey", "UniProt")) %>%
                                arrange(InChIKey, UniProt))))
  
  
  Sd <- Simil.d(list, FP.pub)
  St <- Simil.t(list, target.seq)
  
  #parameters:
  if(dim(Y)[1] < dim(Y)[2]){
    k = as.integer((dim(Y)[1]) / 2)
  } else {
    k = as.integer((dim(Y)[2]) / 2)
  }
  lamb_l = ll
  lamb_d = ld
  lamb_t = lt
  num_iter = 100
  
  #CMF calculates:
  A <- svd(Y)$u[,1:k] %*% (diag(svd(Y)$d)[1:k, 1:k]^0.5)
  B <- svd(Y)$v[,1:k] %*% (diag(svd(Y)$d)[1:k, 1:k]^0.5)
  
  K <- ncol(A)
  lamb_d_Sd <- lamb_d * Sd
  lamb_t_St <- lamb_t * St
  lamb_l_eye_K <- lamb_l * diag(K)
  AtA = t(A) %*% A
  BtB = t(B) %*% B
  
  for(z in 1:num_iter){
    A = (Y %*% B + lamb_d_Sd %*% A) %*% solve((BtB + lamb_l_eye_K + (lamb_d * AtA)))
    AtA = t(A) %*% A
    B = (t(Y) %*% A + lamb_t_St %*% B) %*% solve((AtA + lamb_l_eye_K + (lamb_t * BtB)))
    BtB = t(B) %*% B
  }
  
  y3 <- (A %*% t(B))
  for(p in 1:ncol(y3)){
    y3[ ,p] = normav(y3[ ,p]) 
  }
  
  #reshape the output matrix
  y3 <- reshape2::melt(as.matrix(y3[hclust(dist(y3))$order,
                                    hclust(dist(t(y3)))$order])) %>%
    setNames(c("Drug_ID", "Target_ID", "value")) %>% 
    cbind(method = 'CMF')
  
  return(y3)
  
}

CMF.matrix <- function(list, Sd, St, ll, ld, lt) 
{
  ## list -> dataset input
  ## Sd -> drug similarity matrix
  ## St -> target similarity matrix
  ## ll, ld, lt -> CMF parameters 
  
  #sleecting columns from dataset:
  Y <- as.matrix(table(unique(list %>% select(c("InChIKey", "UniProt")) %>%
                                arrange(InChIKey, UniProt))))
  
  #parameters:
  if(dim(Y)[1] < dim(Y)[2]){
    k = as.integer((dim(Y)[1]) / 2)
  } else {
    k = as.integer((dim(Y)[2]) / 2)
  }
  lamb_l = ll
  lamb_d = ld
  lamb_t = lt
  num_iter = 100
  
  #CMF calculates:
  A <- svd(Y)$u[,1:k] %*% (diag(svd(Y)$d)[1:k, 1:k]^0.5)
  B <- svd(Y)$v[,1:k] %*% (diag(svd(Y)$d)[1:k, 1:k]^0.5)
  
  K <- ncol(A)
  lamb_d_Sd <- lamb_d * Sd
  lamb_t_St <- lamb_t * St
  lamb_l_eye_K <- lamb_l * diag(K)
  AtA = t(A) %*% A
  BtB = t(B) %*% B
  
  for(z in 1:num_iter){
    A = (Y %*% B + lamb_d_Sd %*% A) %*% solve((BtB + lamb_l_eye_K + (lamb_d * AtA)))
    AtA = t(A) %*% A
    B = (t(Y) %*% A + lamb_t_St %*% B) %*% solve((AtA + lamb_l_eye_K + (lamb_t * BtB)))
    BtB = t(B) %*% B
  }
  
  y3 <- (A %*% t(B))
  for(p in 1:ncol(y3)){
    y3[ ,p] = normav(y3[ ,p]) 
  }
  
  #reshape the output matrix
  y3 <- reshape2::melt(as.matrix(y3[hclust(dist(y3))$order,
                                    hclust(dist(t(y3)))$order])) %>%
    setNames(c("Drug_ID", "Target_ID", "value")) %>% 
    cbind(method = 'CMF')
  
  return(y3)
  
}

N_meth <- function(list, list.pred)
{
  
  ## list -> dataset input
  ## list.pred -> dataset predicted
  
  known.inter = unique(list %>% 
                         select(c("InChIKey", "UniProt", "dataset")) %>%
                         setNames(c("Drug_ID", "Target_ID", "dataset")))
  
  list.pred = unique(list.pred %>% 
                       select(c("Target_ID", "Drug_ID", "value", "method", "dataset")))
  
  N_meth = setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                    c("dataset", "method", "cutoff", "Known_DTIs", "New_DTIs"))
  
  AUC_meth = setNames(data.frame(matrix(ncol = 8, nrow = 0)), 
                      c("dataset", "method", "cutoff", "Known_DTIs", "New_DTIs", "AUC_known", "AUC_new", "AUC_rho"))
  
  for (z in unique(list.pred$dataset)) {
    for (x in unique(list.pred$method)){
      for (y in seq(0, 1, 0.05)) {
        
        meth <- unique(list.pred[list.pred$dataset == z & list.pred$method == x,] %>% 
                         filter(value >= y) %>% select(c("Drug_ID", "Target_ID")))
        
        known.meth <- unique(known.inter[known.inter$dataset == z,] %>%
                               select(c("Drug_ID", "Target_ID")))
        
        num.k <- as.numeric((nrow(known.meth) - nrow(anti_join(known.meth, meth))) / nrow(known.meth))
        num.n <- as.numeric((nrow(anti_join(meth, known.meth))) / nrow(meth))
        N_meth <- rbindlist(list(N_meth, as.list(c(z, x, y, num.k, num.n))))
        
      }
      
      N_meth$Known_DTIs <- as.numeric(as.character(N_meth$Known_DTIs))
      N_meth$New_DTIs <- as.numeric(as.character(N_meth$New_DTIs))
      N_meth$cutoff <- as.numeric(as.character(N_meth$cutoff))
      
      v.auc.known <- trapz(N_meth[N_meth$dataset == z & N_meth$method == x,]$cutoff,
                           N_meth[N_meth$dataset == z & N_meth$method == x,]$Known_DTIs)
      v.auc.new <- trapz(N_meth[N_meth$dataset == z & N_meth$method == x,]$cutoff,
                         N_meth[N_meth$dataset == z & N_meth$method == x,]$New_DTIs)
      
      v.auc.rho <- v.auc.new / v.auc.known
      
      temp1 <- N_meth[N_meth$dataset == z & N_meth$method == x,] %>% 
        cbind(AUC_known = v.auc.known, AUC_new = v.auc.new, AUC_rho = v.auc.rho)
      
      AUC_meth <- rbind(AUC_meth, temp1)
      rm(temp1)
      
    }
    
  }
  
  return(AUC_meth)
  
}

corr.target <- function(list, FP.pub, target.seq, met1, met2, ll, ld, lt, r)
{
  
  ## list -> dataset input
  ## FP.pub -> fingerprint from PaDEL(PubChem)
  ## target.seq -> target sequence FASTA (from UniProt)
  ## met1, met2 -> selected methods
  ## ll, ld, lt -> CMF parameters 
  ## r -> rho cutoff
  
  
  # DTI predicted by meth 1
  X.1 <- unique(met1(list, FP.pub, target.seq, ll, ld, lt))
  X.1$Drug_ID <- as.character(X.1$Drug_ID)
  
  # DTI predicted by meth 2
  X.2 <- unique(met2(list, FP.pub, target.seq, ll, ld, lt))
  X.2$Drug_ID <- as.character(X.2$Drug_ID)
  
  corr.target.filter <- setNames(data.frame(matrix(ncol = 6, nrow = 0)),
                                 c("Target_ID", "rho_filter", "rho", "pvalue_spear", 
                                   "cor", "pvalue_pears"))
  
  # loop to calculate spearman correlation by each target
  for (t in unique(list$UniProt)) {
    
    # wraggling the meth1 results only by target t
    target.met1 <- X.1[X.1$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met1", "method"))
    
    # wraggling the meth2 results only by target t
    target.met2 <- X.2[X.2$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met2", "method"))
    
    # suporting the met1 and met2 cluster (by target t)
    # filtering DTI's with zero value in CMF and NBI results
    supp.mets <- cbind(target.met1, target.met2) %>% 
      select(c("Drug_ID", "Target_ID", "value.met1", "value.met2")) %>%
      filter((value.met1 >= 0) | (value.met2 >= 0))
    
    # suporting the met1 and met2 spearman and pearson correlation (by target t)
    supp.spear <- cor.test(supp.mets$value.met1, supp.mets$value.met2, method = "spearman")
    supp.pear <- cor.test(supp.mets$value.met1, supp.mets$value.met2, method = "pearson")
    
    # storing the spearman correlation and p-value by target t
    corr.target.filter <- rbindlist(list(corr.target.filter, 
                                         as.list(c(t, r, supp.spear[[4]][[1]], supp.pear[[4]][[1]], 
                                                   supp.spear[[3]], supp.pear[[3]]))))
    
    rm(supp.spear, supp.pear)
    
  }
  
  # filtering targets with spearman correlation > r
  corr.target.filter <- corr.target.filter %>% 
    filter(rho >= r)
  
  corr.target.filter$rho <- as.numeric(corr.target.filter$rho)
  corr.target.filter$cor <- as.numeric(corr.target.filter$cor)
  
  corr.target.filter$pvalue_spear <- as.numeric(corr.target.filter$pvalue_spear)
  corr.target.filter$pvalue_pears <- as.numeric(corr.target.filter$pvalue_pears)
  
  return(corr.target.filter)
  
}

ensemble.inter <- function(list, FP.pub, target.seq, met1, met2, ll, ld, lt, inter, r)
{
  
  ## list -> dataset input
  ## FP.pub -> fingerprint from PaDEL(PubChem)
  ## target.seq -> target sequence FASTA (from UniProt)
  ## met1, met2 -> selected methods
  ## ll, ld, lt -> CMF parameters 
  ## inter -> top range
  ## r -> rho cutoff
  
  # DTI known
  Y <- unique(list %>% select(c("InChIKey", "UniProt")) %>% 
                setNames(c("Drug_ID", "Target_ID")))
  
  
  # DTI predicted by meth 1
  X.1 <- unique(met1(list, FP.pub, target.seq, ll, ld, lt))
  X.1$Drug_ID <- as.character(X.1$Drug_ID)
  
  # DTI predicted by meth 2
  X.2 <- unique(met2(list, FP.pub, target.seq, ll, ld, lt))
  X.2$Drug_ID <- as.character(X.2$Drug_ID)
  
  # dummy tables to stocking DTIs filtered (target rho >= r) by method
  supp.t1 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                      c("Drug_ID", "Target_ID", "value.met1"))
  
  supp.t2 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                      c("Drug_ID", "Target_ID", "value.met2"))
  
  # loop to calculate spearman correlation by each target
  for (t in unique(Y$Target_ID)) {
    
    # wraggling the meth1 results only by target t
    target.met1 <- X.1[X.1$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met1", "method"))
    
    # wraggling the meth2 results only by target t
    target.met2 <- X.2[X.2$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met2", "method"))
    
    # suporting the met1 and met2 cluster (by target t)
    # filtering DTI's with zero value in CMF and NBI results
    supp.mets <- cbind(target.met1, target.met2) %>% 
      select(c("Drug_ID", "Target_ID", "value.met1", "value.met2")) %>%
      filter((value.met1 >= 0) | (value.met2 >= 0))
    
    # suporting the met1 and met2 spearman correlation (by target t)
    supp.spear <- cor.test(supp.mets$value.met1, supp.mets$value.met2, method = "spearman")
    
    if (supp.spear[[4]][[1]] >= r) {
      
      supp.t1 <- rbind(supp.t1,
                       supp.mets %>% select(c("Drug_ID", "Target_ID", "value.met1")))
      supp.t2 <- rbind(supp.t2,
                       supp.mets %>% select(c("Drug_ID", "Target_ID", "value.met2")))
      
    }
    
    rm(supp.mets)
    
  }
  
  # identify how many new DTIs from meth 1 or meth 2 (filtered)
  
  top = as.numeric(nrow(unique(anti_join(supp.t1, Y, by = c("Drug_ID", "Target_ID")))) + inter)
  
  
  # loop to identify new DTIs intersection from top NBI and top CMF results
  ensem_top <- setNames(data.frame(matrix(ncol = 4, nrow = 0)), 
                        c("Inter_meths", "top", "target_unique", "rho_filter"))
  
  for (e in seq(0, top, inter)) {
    top.inter <- unique(inner_join(unique(anti_join(supp.t1, Y, by = c("Drug_ID", "Target_ID"))) %>%
                                     arrange(-value.met1) %>% dplyr::slice(1:e),
                                   unique(anti_join(supp.t2, Y, by = c("Drug_ID", "Target_ID"))) %>%
                                     arrange(-value.met2) %>% dplyr::slice(1:e),
                                   by = c("Drug_ID", "Target_ID")))
    
    ensem_top <- rbindlist(list(ensem_top, as.list(c(nrow(top.inter), e, length(unique(top.inter$Target_ID)), r))))
    
    rm(top.inter)
    
  }
  
  ensem_top <- ensem_top %>% 
    mutate_if(is.numeric, function(x) ifelse(is.infinite(x), 0, x))
  
  return(ensem_top)
  
  
}

ensemble.DTI <- function(list, FP.pub, target.seq, met1, met2, ll, ld, lt, inter, r, total)
{
  
  ## list -> dataset input
  ## FP.pub -> fingerprint from PaDEL(PubChem)
  ## target.seq -> target sequence FASTA (from UniProt)
  ## met1, met2 -> selected methods
  ## ll, ld, lt -> CMF parameters 
  ## inter -> top range
  ## r -> rho cutoff
  ## total -> total tops
  
  # DTI known
  Y <- unique(list %>% select(c("InChIKey", "UniProt")) %>% 
                setNames(c("Drug_ID", "Target_ID")))
  
  # DTI predicted by meth 1
  X.1 <- unique(met1(list, FP.pub, target.seq, ll, ld, lt))
  X.1$Drug_ID <- as.character(X.1$Drug_ID)
  
  # DTI predicted by meth 2
  X.2 <- unique(met2(list, FP.pub, target.seq, ll, ld, lt))
  X.2$Drug_ID <- as.character(X.2$Drug_ID)
  
  # dummy tables to stocking DTIs filtered (target rho >= r) by method
  supp.t1 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                      c("Drug_ID", "Target_ID", "value.met1"))
  
  supp.t2 <- setNames(data.frame(matrix(ncol = 3, nrow = 0)),
                      c("Drug_ID", "Target_ID", "value.met2"))
  
  # loop to calculate spearman correlation by each target
  for (t in unique(Y$Target_ID)) {
    
    # wraggling the meth1 results only by target t
    target.met1 <- X.1[X.1$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met1", "method"))
    
    # wraggling the meth2 results only by target t
    target.met2 <- X.2[X.2$Target_ID == t, ] %>% arrange(Drug_ID) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met2", "method"))
    
    # suporting the met1 and met2 cluster (by target t)
    # filtering DTI's with zero value in CMF and NBI results
    supp.mets <- cbind(target.met1, target.met2) %>% 
      select(c("Drug_ID", "Target_ID", "value.met1", "value.met2")) %>%
      filter((value.met1 >= 0) | (value.met2 >= 0))
    
    # suporting the met1 and met2 spearman correlation (by target t)
    supp.spear <- cor.test(supp.mets$value.met1, supp.mets$value.met2, method = "spearman")
    
    if (supp.spear[[4]][[1]] >= r) {
      
      supp.t1 <- rbind(supp.t1,
                       supp.mets %>% select(c("Drug_ID", "Target_ID", "value.met1")))
      supp.t2 <- rbind(supp.t2,
                       supp.mets %>% select(c("Drug_ID", "Target_ID", "value.met2")))
      
    }
    
    rm(supp.mets)
    
  }
  
  
  # loop to identify WHO IS new DTIs intersection from top NBI and top CMF results
  ensem_DTIs <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), 
                         c("Drug_ID", "Target_ID", "value.met1", "value.met2", "top"))
  
  
  for (e in seq(10, total, inter)) {
    top.inter <-  unique(inner_join(unique(anti_join(supp.t1, Y, by = c("Drug_ID", "Target_ID"))) %>%
                                      arrange(-value.met1) %>% dplyr::slice(1:e),
                                    unique(anti_join(supp.t2, Y, by = c("Drug_ID", "Target_ID"))) %>%
                                      arrange(-value.met2) %>% dplyr::slice(1:e),
                                    by = c("Drug_ID", "Target_ID"))) 
    
    
    top.inter <- cbind(top.inter, as.list(e)) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met1", "value.met2", "top"))
    
    ensem_DTIs <- rbind(ensem_DTIs, top.inter) %>%
      setNames(c("Drug_ID", "Target_ID", "value.met1", "value.met2", "top"))
    
    rm(top.inter)
    
  }
  
  return(ensem_DTIs)
  
}


######################################################################
############ Just do it! ############ 

###### obtain NBI matrix
NBI.h1 <- NBI(nubbe.ref, FP.pub, target.seq, 2^-2, 2^5, 2^5)


###### obtain drug similarity matrix (Sd)
Sd.halluc <- Simil.d(nubbe.ref, FP.pub)
#Sd.halluc <- Simil.d(halluc.ref, cbind(halluc.ref.drug, halluc.FP.pub) %>% select(-Name))


###### obtain target similarity matrix (St)
St.halluc <- Simil.t(nubbe.ref, target.seq)


###### obtain CMF matrix for the first time
#CMF.h1 <- CMF.matrix(halluc.ref, Sd.halluc, St.halluc, 2^-2, 2^5, 2^5)
CMF.p1 <- CMF(nubbe.ref, FP.pub, target.seq, 2^-2, 2^5, 2^5)


###### sensibility method - NBI
AUC.p1.NBI <- N_meth(nubbe.ref %>% cbind(dataset = 'Plants'),
                     NBI.p1 %>% cbind(dataset = 'Plants'))


###### sensibility method - CMF
AUC.p1.CMF <- N_meth(nubbe.ref %>% cbind(dataset = 'Plants'),
                     CMF.h1 %>% cbind(dataset = 'Plants'))


###### spearman correlation : NBI/CMF (filtering)
corr.spear.p1 <- corr.target(nubbe.ref, FP.pub, target.seq, NBI, CMF, 2^-2, 2^5, 2^5, 0.5) %>%
  cbind(dataset = 'Plants')


###### spearman correlation : NBI/CMF (NO filter)
corr.spear.nofilter.p1 <- corr.target(nubbe.ref, FP.pub, target.seq, NBI, CMF, 2^-2, 2^5, 2^5, 0) %>%
  cbind(dataset = 'Plants (no filter)')


###### NBI and CMF common DTIs x tops (scatter) 
ensemble.inter.p1 <- ensemble.inter(nubbe.ref, FP.pub, target.seq, NBI, CMF, 2^-2, 2^5, 2^5, 10, 0.5) %>%
  cbind(dataset = 'Plants')


###### NBI and CMF common DTIs x frequency (table) 
ensemble.DTI.p1 <- ensemble.DTI(nubbe.ref, FP.pub, target.seq, NBI, CMF, 2^-2, 2^5, 2^5, 10, 0.5, 50) %>%
  cbind(dataset = 'Plants')








