#' This function merge the flyID by either genesymbol, CG_ID or flybase ID
#'
#' @param input data.frame to merge with 
#' @param merge_id_x colname in input data frame to merge with 
#' @param merge_id_y colname in the reference file, one of "CG_ID",  "flybaseID" , "genesymbol"
#' @export
#' @import dplyr 
#' @import tidyr
#' @examples
#' \dontrun{
#' input = structure(list(id1 = c("FBgn0011656", "FBgn0011656", "FBgn0011656", "FBgn0011656", "FBgn0011656", "FBgn0011656"), id2 = c("FBgn0004009", "FBgn0004110", "FBgn0010246", "FBgn0039039", "FBgn0086906", "FBgn0264491")), .Names = c("id1", "id2"), row.names = c(NA, -6L), class = "data.frame")
#' input%>%mergeIDs(input=.,merge_id_x="id1",merge_id_y="flybaseID")%>%mergeIDs(input=.,merge_id_x="id2",merge_id_y="flybaseID")%>%arrange(id2,type)%>%filter(type!="rna")
#' }


mergeIDs= function(input,merge_id_x = c('id1'),merge_id_y=c("CG_ID",  "flybaseID" , "genesymbol"))
{
data(FlyBase_IDs)
ID_data=  FlyBase_IDs%>%dplyr::rename(CG_ID=submitted_id,flybaseID = current_id, genesymbol = current_symbol)%>%select(-converted_id )%>%unique
names(ID_data) = paste(merge_id_x,names(ID_data),sep="_")

output = input%>%merge(.,ID_data,by.x=merge_id_x,by.y=paste(merge_id_x,merge_id_y,sep="_"),all.x=TRUE)

return(output)
}



