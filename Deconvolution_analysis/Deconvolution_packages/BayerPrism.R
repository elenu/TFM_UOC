# https://github.com/Danko-Lab/BayesPrism

library("devtools");
install_github("Danko-Lab/BayesPrism/BayesPrism")

suppressWarnings(library(BayesPrism))
setwd("~/Desktop/Elena_R_scripts/Deconvolution")

# Load data
load("tutorial.gbm.rdata")
ls()
dim(bk.dat)
head(rownames(bk.dat))
head(colnames(bk.dat))
dim(sc.dat)
head(rownames(sc.dat))
head(colnames(sc.dat))
sort(table(cell.type.labels))
sort(table(cell.state.labels))
table(cbind.data.frame(cell.state.labels, cell.type.labels))

# QC of data
plot.cor.phi (input=sc.dat,
              input.labels=cell.state.labels,
              title="cell state correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.cs",
              cexRow=0.2, cexCol=0.2,
              margins=c(2,2))
plot.cor.phi (input=sc.dat,
              input.labels=cell.type.labels,
              title="cell type correlation",
              #specify pdf.prefix if need to output to pdf
              #pdf.prefix="gbm.cor.ct",
              cexRow=0.5, cexCol=0.5,
)

# Filter outlier genes
sc.stat <- plot.scRNA.outlier(
  input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE #return the data used for plotting.
  #pdf.prefix="gbm.sc.stat" specify pdf.prefix if need to output to pdf
)
head(sc.stat) 
bk.stat <- plot.bulk.outlier(
  bulk.input=bk.dat,#make sure the colnames are gene symbol or ENSMEBL ID
  sc.input=sc.dat, #make sure the colnames are gene symbol or ENSMEBL ID
  cell.type.labels=cell.type.labels,
  species="hs", #currently only human(hs) and mouse(mm) annotations are supported
  return.raw=TRUE
  #pdf.prefix="gbm.bk.stat" specify pdf.prefix if need to output to pdf
)
head(bk.stat)
sc.dat.filtered <- cleanup.genes (input=sc.dat,
                                  input.type="count.matrix",
                                  species="hs",
                                  gene.group=c( "Rb","Mrp","other_Rb","chrM","MALAT1","chrX","chrY") ,
                                  exp.cells=5)
dim(sc.dat.filtered)
plot.bulk.vs.sc (sc.input = sc.dat.filtered,
                 bulk.input = bk.dat
                 #pdf.prefix="gbm.bk.vs.sc" specify pdf.prefix if need to output to pdf
)
sc.dat.filtered.pc <- select.gene.type (sc.dat.filtered,
                                        gene.type = "protein_coding")
# Select marker genes (Optional)
# performing pair-wise t test for cell states from different cell types
diff.exp.stat <- get.exp.stat(sc.dat=sc.dat[,colSums(sc.dat>0)>3],# filter genes to reduce memory use
                              cell.type.labels=cell.type.labels,
                              cell.state.labels=cell.state.labels,
                              pseudo.count=0.1, #a numeric value used for log2 transformation. =0.1 for 10x data, =10 for scell.count.cutoff=50, # a numeric value to exclude cell state with number of cells fewer thann.cores=1 #number of threads
)
sc.dat.filtered.pc.sig <- select.marker (sc.dat=sc.dat.filtered.pc,
                                         stat=diff.exp.stat,
                                         pval.max=0.01,
                                         lfc.min=0.1)
dim(sc.dat.filtered.pc.sig)

# Contruct a Prism object
myPrism <- new.prism(
  reference=sc.dat.filtered.pc,
  mixture=bk.dat,
  input.type="count.matrix",
  cell.type.labels = cell.type.labels,
  cell.state.labels = cell.state.labels,
  key="tumor",
  outlier.cut=0.01,
  outlier.fraction=0.1,
)

# Run Bayes Prim
bp.res <- run.prism(prism = myPrism, n.cores=3) ## Messages include: "estimated time to complete: 45mins"
bp.res

# Extract results
slotNames(bp.res)
# extract posterior mean of cell type fraction theta
theta <- get.fraction (bp=bp.res,
                       which.theta="final",
                       state.or.type="type")
head(theta)
# extract coefficient of variation (CV) of cell type fraction
theta.cv <- bp.res@posterior.theta_f@theta.cv
head(theta.cv)
# extract posterior mean of cell type-specific gene expression count matrix Z
Z.tumor <- get.exp (bp=bp.res,
                    state.or.type="type",
                    cell.name="tumor")
head(t(Z.tumor[1:5,]))
#save the result
save(bp.res, file="bp.res.rdata")
