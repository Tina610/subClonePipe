arg<-commandArgs(T)

if(length(arg) != 4){
  print('Rscript this.R cloneEvaInput.txt Outdir sump alpha')
  q()
}

library(clonevol)
library(fishplot)
createFishPlotObjects <- function(results){
  fishes = list()
  for(i in 1:results$num.models){
    print(i)
    fishes[[i]]=createFishObject(results$cell.fractions[[i]],
                                 parents=as.integer(results$parents[[i]]),
                                 timepoints=results$timepoints,
                                 clone.labels=results$clonevol.clone.names[1:length(results$clonevol.clone.names)],
                                 fix.missing.clones=TRUE)
  }
  return(fishes)
}
# read data
data<-read.table(arg[1], header=TRUE,sep='\t', quote='')
# sort data
data <- data[order(data$cluster),]
#########################################
ncol = length(names(data))
usenames = names(data)[-c(1:3,ncol)]
n = length(usenames)/2
vaf.col.names=usenames[1:n]
sample.groups <- vaf.col.names;
#########################################
names(sample.groups) <- vaf.col.names
data[vaf.col.names]=data[vaf.col.names]*100
y = infer.clonal.models(variants = data,
                        cluster.col.name = 'cluster',
                        vaf.col.names = vaf.col.names,
                        sample.groups = sample.groups,
                        cancer.initiation.model='monoclonal',
                        subclonal.test = 'bootstrap',
                        subclonal.test.model = 'non-parametric',
                        num.boots = 1000,
                        founding.cluster = NULL,
                        cluster.center = 'mean',
                        ignore.clusters = NULL,
                        #clone.colors = clone.colors,
                        min.cluster.vaf = 0,
                        # min probability that CCF(clone) is non-negative
                        sum.p = args[3],
                        # alpha level in configendence interval estimate for CCF(clone)
                        alpha = args[4])
print(y)

cloeva = y

y <- transfer.events.to.consensus.trees(y,
                                        data[data$Driver,],
                                        cluster.col.name = 'cluster',
                                        event.col.name = 'id')

y <- convert.consensus.tree.clone.to.branch(y, branch.scale = 'sqrt')

plot.clonal.models(y,
                   # box plot parameters
                   box.plot = TRUE,
                   fancy.boxplot = TRUE,
                   fancy.variant.boxplot.vaf.limits = 100,
                   fancy.variant.boxplot.highlight = 'Driver',
                   fancy.variant.boxplot.highlight.shape = 21,
                   fancy.variant.boxplot.highlight.fill.color = 'red',
                   fancy.variant.boxplot.highlight.color = 'black',
                   fancy.variant.boxplot.highlight.note.col.name = 'id',
                   fancy.variant.boxplot.highlight.note.color = 'blue',
                   fancy.variant.boxplot.highlight.note.size = 2,
                   fancy.variant.boxplot.jitter.alpha = 1,
                   fancy.variant.boxplot.jitter.center.color = 'grey50',
                   fancy.variant.boxplot.base_size = 12,
                   fancy.variant.boxplot.plot.margin = 1,
                   fancy.variant.boxplot.vaf.suffix = '.VAF',
                   # bell plot parameters
                   clone.shape = 'bell',
                   bell.event = TRUE,
                   bell.event.label.color = 'blue',
                   bell.event.label.angle = 60,
                   clone.time.step.scale = 1,
                   bell.curve.step = 2,
                   # node-based consensus tree parameters
                   merged.tree.plot = TRUE,
                   tree.node.label.split.character = NULL,
                   tree.node.shape = 'circle',
                   tree.node.size = 30,
                   tree.node.text.size = 0.5,
                   merged.tree.node.size.scale = 1.25,
                   merged.tree.node.text.size.scale = 2.5,
                   merged.tree.cell.frac.ci = FALSE,
                   # branch-based consensus tree parameters
                   merged.tree.clone.as.branch = TRUE,
                   mtcab.event.sep.char = ',',
                   mtcab.branch.text.size = 1,
                   mtcab.branch.width = 0.75,
                   mtcab.node.size = 3,
                   mtcab.node.label.size = 1,
                   mtcab.node.text.size = 1.5,
                   # cellular population parameters
                   cell.plot = TRUE,
                   num.cells = 100,
                   cell.border.size = 0.25,
                   cell.border.color = 'black',
                   clone.grouping = 'horizontal',
                   #meta-parameters
                   scale.monoclonal.cell.frac = TRUE,
                   show.score = FALSE,
                   cell.frac.ci = TRUE,
                   disable.cell.frac = FALSE,
                   # output figure parameters
                   out.dir = arg[2],
                   out.format = 'pdf',
                   overwrite.output = TRUE,
                   width = 8,
                   height = 4,
                   # vector of width scales for each panel from left to right
                   panel.widths = c(3,4,2,4,2))
#dev.off()


print('fishplot')

f = generateFishplotInputs(results=cloeva)
fishes = createFishPlotObjects(f)
#plot with fishplot
pdf(paste(arg[2],'/fish.pdf',sep=''), width=8, height=5)
for (i in 1:length(fishes)){
  fish = layoutClones(fishes[[i]])
  fish = setCol(fish,f$clonevol.clone.colors)
  fishPlot(fish,shape="spline", title.btm="Patient", cex.title=0.5,
           vlines=seq(1, length(sample.groups)), vlab=sample.groups, pad.left=0.5)
}
dev <- dev.off()





