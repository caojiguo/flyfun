#' This function identity search the interactions in the flybase data
#'
#' @param input_interactions data.frame with two cols
#' @param s interaction type one of "ppi","rna","tf","generic"
#' @import dplyr 
#' @import tidyr

search_each = function(s,input_interactions=input_int)
{

names(input_interactions) = c('id1','id2')	

if (s=="ppi")
{
data(flybase_ppi)
flybase_ppi = flybase_ppi%>%dplyr::rename(id1 = FLY_GENE1,id2=FLY_GENE2)
input_interactions2 = input_interactions[,c("id2","id1")];
names(input_interactions2) = c('id1','id2')
input_interactions = rbind.data.frame(input_interactions,input_interactions2)
res = input_interactions%>%merge(.,flybase_ppi,by=c("id1","id2"))
}

if (s=="rna")
{

data(rna_gene)
fun <- function(x, y){ 
rna1 = as.character(subset(rna_gene,FLY_TARGET_GENE%in%x)$RNA_SYMBOL)
rna2 = as.character(subset(rna_gene,FLY_TARGET_GENE%in%y)$RNA_SYMBOL)
data.frame(rna=intersect(rna1,rna2))
}

res = input_interactions%>%rowwise()%>%do(rna = fun(.$id1[1],.$id2[1]))%>%mutate(n=nrow(rna),rna = paste0(rna[,1],collapse="\\"))%>%select(n,rna)%>%cbind.data.frame(input_interactions,.)


}

if (s=="tf")
{
data(tf_gene)
tf_gene= tf_gene%>%dplyr::rename(id1 = FLY_TF_GENE,id2=FLY_TARGET_GENE)
res = input_interactions%>%merge(.,tf_gene,by=c("id1","id2"))
}

if (s=="genetic")
{

data(fly_genetic_interactions)
fly_genetic_interactions= fly_genetic_interactions%>%rename(id1 = FLY_GENE1,id2=FLY_GENE2)
input_interactions2 = input_interactions[,c("id2","id1")];
names(input_interactions2) = c('id1','id2')
input_interactions = rbind.data.frame(input_interactions,input_interactions2)
res = input_interactions%>%merge(.,fly_genetic_interactions,by=c("id1","id2"))
}

return(res)
}
