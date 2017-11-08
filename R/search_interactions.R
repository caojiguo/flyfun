#' This function identity search the interactions in the flybase data
#'
#' @param input_interactions data.frame with two cols
#' @param type interaction type one of 'all',"ppi","rna","tf","generic"
#' @export
#' @import dplyr 
#' @import tidyr
#' @examples
#' \dontrun{
#' input_interactions = data.frame(id1 = c('FBgn0011656','FBgn0003149',"7872"),id2 = c('FBgn0003149','FBgn0003900',"9232"))
#' output = search_interaction(input_interactions=pairs,type="all")
#' output%>%do(data.frame(type=.$type,.$outdata%>%select(id1,id2)))%>%group_by(type)%>%summarise(n=n())
#' output%>%do(data.frame(type=.$type,.$outdata))%>%select(id1,id2,PMID_URL)%>%tail
#' }

search_interaction = function(input_interactions,type= c('all',"ppi","rna","tf","generic"))
{


input_int = unique(input_interactions)

source = type

if (type=="all") { source = c("ppi","rna","tf","genetic") } 

output = source%>%data.frame(type=.)%>%group_by(type)%>%do(outdata=search_each(.$type,input_interactions=input_int))


return(output)

}



