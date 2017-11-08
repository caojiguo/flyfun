#' This function plot the gene trajectory from drosophila data at embryonic stage 
#'
#' @param input_ids a vector contains either flybase ID, CG_ID or genesymbol; can be mixed
#' @export
#' @import dplyr 
#' @import tidyr
#' @import ggplot2
#' @importFrom plyr ldply
#' @examples
#' \dontrun{
#' input_ids = c("CG11880", "tin", "CG6404","XXXXXXX","d4", "Traf6")
#' gene_curve(input_ids)
#' }

gene_curve= function(input_ids,plot.only = TRUE)
{
data(drosophila)
data(timepoints)
res = ldply(input_ids, function(input){

print(input)
CG_ID = try(matched_id(input)$CG_ID[1],silent=TRUE)
if(class(CG_ID)=="try-error") {
 	print(' no matched CG ID; please enter the valiad gene name.')
 	return(NULL)
} else {
	dataplot = data.frame(value=2^drosophila[CG_ID,timepoints$symbol],time=timepoints$realtime)
	dataplot$ID= paste0(matched_id(input)%>%head(.,1)%>%as.character%>%unique,collapse="_")
	
	return(dataplot)
}

}
)
plot  = ggplot(res,aes(x=time,y=value))+geom_line()+geom_smooth(se=FALSE)+facet_wrap(~ID,scales='free')+theme_bw()
if(plot.only) return(plot) else  return(list(plot=plot, data=res))
}





