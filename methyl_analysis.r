library(methylKit)

file.list=list('sample1_CG.txt',
               'sample2_CG.txt',
               'sample3_CG.txt',
               'sample4_CG.txt',
               'sample5_CG.txt'
)

myobj_5 = methRead(location = file.list, sample.id = list("sample1","sample2","sample3", "sample4", "sample5"), 
                   assembly = "hg38", dbtype = "tabix", pipeline = "bismarkCytosineReport",
                   header = FALSE, skip = 0, sep = "\t", resolution = "region", 
                   treatment = c(0,0,1,1,1), dbdir = "Lung_DB", mincov = 5,)

tile500_5 = tileMethylCounts(myobj_5, win.size=500, step.size=500, cov.bases = 0, mc.cores = 1)
tile500_5filtered = filterByCoverage(tile500_5, lo.count = 5, lo.perc = NULL, hi.count = NULL, hi.perc = 99.9)
tile500_5normalized = normalizeCoverage(tile500_5filtered, method = "median")
tileMeth5 = unite(tile500_5normalized)

pm5=percMethylation(tileMeth5)
mds5=matrixStats::rowSds(pm5)

names(mds5)<- 1:length(mds5)

mds5

hist(mds5,col="grey",xlab="Std.dev. per CpG")

pdf(file="LungCorrelationPlot_tile500_5.pdf",
    width=12, height=12); getCorrelation(tileMeth5, plot=TRUE); dev.off()

pdf(file="Lung.ctreeplot_tile500_5.pdf", 
    width=12, height=12); clusterSamples(tileMeth5, dist="correlation", method="ward", plot=TRUE); dev.off()

pdf(file="Lung.screeplot_tile500_5.pdf",
    width=12,height=12); PCASamples(tileMeth5, screeplot=TRUE); dev.off()

pdf(file="Lung.pcaxyplot_tile500_5.pdf",
    width=12,height=12); PCASamples(tileMeth5); dev.off()

myDiff_tile500_5=calculateDiffMeth(tileMeth5)

myDiff25p_hyper_tile500_5 <- getMethylDiff(myDiff_tile500_5, difference=25, qvalue=0.01, type= "hyper",)
myDiff25p_hypo_tile500_5 <- getMethylDiff(myDiff_tile500_5, difference=25, qvalue=0.01, type="hypo",)

#pdf("DiffMethPerChr_tile500.pdf",width=12,height=12); 
#diffMethPerChr(myDiff_tile500_10,plot=TRUE,qvalue.cutoff=0.01,meth.cutoff=25,
# exclude = c("chr6_ssto_hap7","chr6_mcf_hap5","chr6_cox_hap2","chr6_mann_hap4","chr6_apd_hap1",
#"chr6_qbl_hap6","chr6_dbb_hap3","chr17_ctg5_hap1","chr4_ctg9_hap1","chr1_gl000192_random",
# "chrUn_gl000225","chr4_gl000194_random","chr4_gl000193_random","chr9_gl000200_random",
# "chrUn_gl000222","chrUn_gl000212","chr7_gl000195_random","chrUn_gl000223",
# "chrUn_gl000224","chrUn_gl000219","chr17_gl000205_random","chrUn_gl000215",
# "chrUn_gl000216","chrUn_gl000217","chr9_gl000199_random","chrUn_gl000211",
# "chrUn_gl000213","chrUn_gl000220","chrUn_gl000218","chr19_gl000209_random",
# "chrUn_gl000221","chrUn_gl000214","chrUn_gl000228","chrUn_gl000227",
# "chr1_gl000191_random","chr19_gl000208_random","chr9_gl000198_random",
# "chr17_gl000204_random","chrUn_gl000233","chrUn_gl000237","chrUn_gl000230",
# "chrUn_gl000242","chrUn_gl000243","chrUn_gl000241","chrUn_gl000236",
# "chrUn_gl000240","chr17_gl000206_random","chrUn_gl000232","chrUn_gl000234",
# "chr11_gl000202_random","chrUn_gl000238","chrUn_gl000244","chrUn_gl000248",
# "chr8_gl000196_random","chrUn_gl000249","chrUn_gl000246","chr17_gl000203_random",
# "chr8_gl000197_random","chrUn_gl000245","chrUn_gl000247","chr9_gl000201_random",
# "chrUn_gl000235","chrUn_gl000239","chr21_gl000210_random","chrUn_gl000231",
# "chrUn_gl000229","chrM","chrUn_gl000226","chr18_gl000207_random")); dev.off()


bedgraph(myDiff25p_hyper_tile500_5,file.name="Methdiff.bedgraph_hyper_tile500_5",col.name="meth.diff")
bedgraph(myDiff25p_hypo_tile500_5,file.name="Methdiff.bedgraph_hypo_tile500_5",col.name="meth.diff")

library(genomation)

gene.obj <- genomation::readTranscriptFeatures("UCSC_hg38_genepredbed12.bed")

diffCpG_regions_ann_hypo_5 = annotateWithGeneParts(as(myDiff25p_hypo_tile500_5,"GRanges"), gene.obj)
diffCpG_regions_ann_hyper_5 = annotateWithGeneParts(as(myDiff25p_hyper_tile500_5,"GRanges"), gene.obj)

pdf("hypoCpG_annotate_tile500_5.pdf",
    width=12,height=12); plotTargetAnnotation(diffCpG_regions_ann_hypo_5); dev.off()
pdf("hyperCpG_annotate_tile500_5.pdf",
    width=12,height=12); plotTargetAnnotation(diffCpG_regions_ann_hyper_5); dev.off()

write.table(getAssociationWithTSS(diffCpG_regions_ann_hypo_5), file="TSS_hypo_tile500_5.txt")
write.table(getAssociationWithTSS(diffCpG_regions_ann_hyper_5), file="TSS_hyper_tile500_5.txt")

tss.assoc_hypo <- getAssociationWithTSS(diffCpG_regions_ann_hypo_5)
tss.assoc_hyper <- getAssociationWithTSS(diffCpG_regions_ann_hyper_5)

pdf("TSS_dist_hypoCpG_tile500_5.pdf", width=12,height=12); hist(tss.assoc_hypo$dist.to.feature[abs(tss.assoc_hypo$dist.to.feature)<=100000],
                                                                main="Distance to nearest TSS", xlab="Distance in basepairs", breaks=50, col="brown4"); dev.off()

pdf("TSS_dist_hyperCpG_tile500_5.pdf",width=12,height=12); hist(tss.assoc_hyper$dist.to.feature[abs(tss.assoc_hyper$dist.to.feature)<=100000],
                                                                main="Distance to nearest TSS", xlab="Distance in basepairs", breaks=50, col="brown4"); dev.off()

write.table(getMembers(diffCpG_regions_ann_hyper_5), file="CpG_hypermeth_tile500_5.txt")
write.table(getMembers(diffCpG_regions_ann_hypo_5), file="CpG_hypometh_tile500_5.txt")

CpG_region_hypo <- tss.assoc_hypo[,1]
CpG_region_hyper <- tss.assoc_hyper[,1]

write.table(myDiff25p_hypo_tile500_5[CpG_region_hypo,], file="Hypo_cpg_region_tile500_5.tsv", sep='\t', quote=FALSE)
write.table(myDiff25p_hyper_tile500_5[CpG_region_hyper,], file="Hyper_cpg_region_tile500_5.tsv", sep='\t', quote=FALSE)

#annotation of differentially methylated hyper and hypo CpG
cpg.obj <- readFeatureFlank("reg_CpG_hg38.bed.txt", feature.flank.name=c("CpGi","Shores"))

diffCpGann_hyper_tile500 = annotateWithFeatureFlank(as(myDiff25p_hyper_tile500_5,"GRanges"),
                                                    cpg.obj$CpGi,cpg.obj$Shores,
                                                    feature.name="CpGi",flank.name="Shores")

write.table(getMembers(diffCpGann_hyper_tile500), file="hyper_CpGi_shore_5.txt")
plotTargetAnnotation(diffCpGann_hyper_tile500, main = "Hypermethylated-CpG annotation")

diffCpGann_hypo_tile500 = annotateWithFeatureFlank(as(myDiff25p_hypo_tile500_5,"GRanges"),
                                                   cpg.obj$CpGi,cpg.obj$Shores,
                                                   feature.name="CpGi",flank.name="Shores")

write.table(getMembers(diffCpGann_hypo_tile500), file="hypo_CpGi_shore_5.txt")
plotTargetAnnotation(diffCpGann_hypo_tile500, main = "Hypomethylated-CpG annotation")

CpG_annotation_hyper <- read.table("hyper_CpGi_shore_5.txt", header=T, stringsAsFactors = F, row.names = 1)
CpG_annotation_hypo <- read.table("hypo_CpGi_shore_5.txt", header=T, stringsAsFactors = F, row.names = 1)

#adding exon promoter and intron info in the annotation for hypermethylated CpG
hypermeth_CpG <- getData(myDiff25p_hyper_tile500_5)

annotateHyperCpG <- merge(hypermeth_CpG, CpG_annotation_hyper, by=0, all=TRUE)

hyperdiff_exon_prom_intron <- read.table("CpG_hypermeth_tile500_5.txt", header=T, stringsAsFactors = F, row.names = 1)

hyper_annotate <- cbind(annotateHyperCpG, hyperdiff_exon_prom_intron)
write.table(hyper_annotate, file="hyper_annotated_CpG_5.txt",sep='\t', quote=FALSE)


#adding exon promoter and intron info in the annotation for hypomethylated CpG
hypometh_CpG <- getData(myDiff25p_hypo_tile500_5)

annotateHypoCpG <- merge(hypometh_CpG, CpG_annotation_hypo, by=0, all=TRUE)

hypodiff_exon_prom_intron <- read.table("CpG_hypometh_tile500_5.txt", header=T, stringsAsFactors = F, row.names = 1)

hypo_annotate <- cbind(annotateHypoCpG, hypodiff_exon_prom_intron)
write.table(hypo_annotate, file="hypo_annotated_CpG_5.txt",sep='\t', quote=FALSE)


##############################################################################################################################################
##############################################################################################################################################
# Heatmap for hypomethylated regions
# here, 1:1160 range comes from the object's size; vary according to your object size.
query_hypo_rows<-paste(myDiff25p_hypo_tile500_5[1:1160]$chr,":",myDiff25p_hypo_tile500_5[1:1160]$start, "-", myDiff25p_hypo_tile500_5[1:1160]$end, sep = "")
#dataframe for all DMLoci
all_loci_samples <- getData(tileMeth5)
# placing row name as chr:start
rownames(all_loci_samples)<- paste(all_loci_samples$chr, ":", all_loci_samples$start, "-", all_loci_samples$end, sep = "")
# getting only those lines from all_loci_samples that are present in query_hyper_rows
query_hypo_res<- all_loci_samples[rownames(all_loci_samples) %in% query_hypo_rows,]
#made a copy
query_hypo_res_natozero <- query_hypo_res
#replacing na with 0
#query_hypo_res_natozero[is.na(query_hypo_res_natozero)] <- 0
#dividing col having Cs to Coverage and storing in a columns meth1,2,3...so on
#change the below columns as per your needs
query_hypo_res_natozero$meth1 <- query_hypo_res[,6]/query_hypo_res[,5]
query_hypo_res_natozero$meth2 <- query_hypo_res[,9]/query_hypo_res[,8]
query_hypo_res_natozero$meth3 <- query_hypo_res[,12]/query_hypo_res[,11]
query_hypo_res_natozero$meth4 <- query_hypo_res[,15]/query_hypo_res[,14]
query_hypo_res_natozero$meth5 <- query_hypo_res[,18]/query_hypo_res[,17]


#again converted na to zero because above computation generated NA values
#query_hypo_res_natozero[is.na(query_hypo_res_natozero)] <- 0
#round the ratio to 2 decimal place; it just rewritten col 29 to 36 in the dataframe while removing other columns
query_hypo_res_natozero<-round(query_hypo_res_natozero[,(16:20)], digits = 2)
#renaming the columns in the dataframe just created above
colnames(query_hypo_res_natozero)[1] = "sample1"
colnames(query_hypo_res_natozero)[2] = "sample2"
colnames(query_hypo_res_natozero)[3] = "sample3"
colnames(query_hypo_res_natozero)[4] = "sample4"
colnames(query_hypo_res_natozero)[5] = "sample5"


#saves the dataframe to a file

write.table(query_hypo_res_natozero, file="query_hypo_res_natozero_5.tsv", sep='\t', quote=FALSE)

#Annotation_columns creation
condition<- c("Lung cancer","Lung cancer","Lung cancer","Lung cancer","Lung cancer","Lung cancer",
              "Control","Control","Control","Control","Control", "Control")
sex <- c("female","male","male", "male", "male", "male", "na","female","female","male","male","male")

#created a dataframe patients
patients <- data.frame(condition,sex)
data.frame(patients, stringsAsFactors = True)

rownames(patients)<- c("BTD083","ACRI500","BTD075", "BTD097", "BTD122", "BTD159", 
                       "ACRI002","ACRI069", "ACRI075", "Prost015", "Prost018", "Prost022")
write.table(patients, file="ALSpatients.csv", sep=',', quote=FALSE)

library(pheatmap)
#library(heatmaply)
CpGloci_hypo <- read.csv("query_hypo_res_natozero_5.tsv", header=T, sep='\t',
                         row.names = 1)

#annotation file for column
annotation_col <- read.csv("ALSpatients.csv", header = T, row.names = 1)

#creating annotation_row file
dt<- read.table("CpG_hypometh_tile500_5.txt", header = T, stringsAsFactors = F, row.names = 1)
head(dt)
v <- ifelse(dt$prom == 1, "Promoter", "Other")
v[dt$exon == 1]<- "Exon"
v[dt$intron == 1]<- "Intron"
dt$new <- v
head(dt)
#write.table(dt, file="diffCpGregion_hypo.csv", sep=',', quote=FALSE)
#"diffCpGregion_hyper_v6.csv" has more column and i only wanted "new" column 
#with Promoter/exon/Intron/Others, so i did

copy_dt <- dt

#loading dyplr as there was an error for %>% not found and this package will load it
library(dplyr)

#just adding the correcsponding position chr:start as row name in my copy_dt.
#did not shuffle the dataframe, just pasted the sequence of chr:start in the
#same order as it is present in the query_hyper_res_natozero
rownames(copy_dt)<-paste(rownames(query_hypo_res_natozero))

#similarly copied the action in dt
rownames(dt)<-paste(rownames(query_hypo_res_natozero))

#changed the rowname from "new" to "region"
colnames(copy_dt)[4] = "region"
undesired <- c("prom","exon","intron")
copy_dt_regions <- copy_dt %>% dplyr::select(-one_of(undesired))
write.table(copy_dt_regions, file="diffCpGregion_hypo_5.csv", sep=',', quote=FALSE)

#Annoation file for row
annotation_row <- read.csv("diffCpGregion_hypo_5.csv", header = T,
                           row.names = 1)

str(annotation_row)
str(CpGloci_hypo)
str(annotation_col)
annotation_col

#color coding for the legends
annotation_colors <- list(condition = c('Lung cancer'="cornflowerblue", 'Control'= "#FC9272"),
                          sex = c('female' = "orange", 'male' = "brown", 'na' = "grey"), 
                          region = c('Promoter'="blue", 'Exon'="purple", 'Intron'="pink",
                                     'Other' = "light blue"))

pdf("heatmap_hypoCpG_tile500_5.pdf",width=12,height=12); pheatmap(CpGloci_hypo, fontsize = 12, fontsize_row = 0.1, fontsize_col = 8, cluster_rows = T, 
                                                                  cluster_cols = T, annotation_row = annotation_row, annotation_col = annotation_col, 
                                                                  cellwidth = 25, cellheight = 0.5, angle_row = 45,
                                                                  main = "Hypomethylation across Lung Cancer and Control samples", angle_col = 45,
                                                                  treeheight_row = 100, treeheight_col = 50, annotation_colors = annotation_colors, na_col = "grey",
                                                                  color = colorRampPalette(c("white", "yellow", "green", "blue"))(100)); dev.off()




##########Similarly do for Heatmap for hypermethylated regions###########

