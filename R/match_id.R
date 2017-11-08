#' This function identity the flybase ID, CG_ID or genesymbol
#'
#' @param input_ids a vector contains either flybase ID, CG_ID or genesymbol; can be mixed
#' @export
#' @import dplyr 
#' @import tidyr
#' @importFrom plyr ldply
#' @examples
#' \dontrun{
#' input_ids = c("CG11880", "tin", "CG6404","XXXXXXX","d4", "Traf6")
#' matched_id(input_ids)
#' }

matched_id = function(input_ids)
{
data(FlyBase_IDs)
ID_data=  FlyBase_IDs%>%dplyr::rename(CG_ID=submitted_id,flybaseID = current_id, genesymbol = current_symbol)%>%dplyr::select(-converted_id)%>%unique
allIDs= ID_data%>%unite(allnames,CG_ID,flybaseID,genesymbol,sep=";")
output0 =  plyr::ldply(input_ids,function(x) data.frame(input = x ,match = any(grepl(x,allIDs$allnames,fixed=TRUE))))


output = plyr::ldply(output0$input[output0$match],function(x) data.frame(input = x ,match = grep(x,allIDs$allnames,value=TRUE,
fixed=TRUE)))

output%>%separate(match,c("CG_ID","flybaseID","genesymbol"),sep=";")%>%dplyr::filter(input==CG_ID|input==flybaseID|input==genesymbol)


}


