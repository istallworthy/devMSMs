#' Wide complete data (continuous exposure)
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). These data contain economic strain (ESEATA1) as a continuously 
#' distributed variable. 
#'
#' @name sim_data_wide.rda
#' @docType data
#' @format A wide data frame of 1,292 observations
#' There are 69 measured variables. 
#' \itemize{
#' \item "ESETA1" is the continuous exposure of economic strain
#' \item "StrDif_Tot.58" is the continuous outcome of behavioral problems
#' \item "InRatioCor" is the income-to-needs ratio
#' \item "PmEd2" is the parent's education level
#' \item "state" is the family's state of residence
#' \item "TcBlac2" is the family's race (1 = x, 0 = y)
#' \item "bioDadInHH2" is whether the biological father lives with the family (insert coding)
#' \item "HomeOwnd" indicator of whether family owns home (insert coding)
#' \item "KFASTScr"
#' \item "PmBlac2" primary careigver race (insert coding)
#' \item "SmokTotl"
#' \item "caregiv_health"
#' \item "gov_assist"
#' \item "ALI_LE"
#' \item "B18Raw"
#' \item "CORTB"
#' \item "EARS_TJo"
#' \item "fscore"
#' \item "HOMEETA1"
#' \item "IBRAttn"
#' \item "LESMnNeg"
#' \item "MDI"
#' \item "RHAsSO"
#' \item "SAAmylase"
#' \item "WndNbrhood"
#' }
#' @references Vernon-Feagans, L., Cox, M., Willoughby, M., Burchinal, M., Garrett-Peters, P., Mills-Koonce, R., 
#' Garrett-Peiers, P., Conger, R. D., & Bauer, P. J. (2013). The Family Life Project: An Epidemiological and 
#' Developmental Study of Young Children Living in Poor Rural Communities.
#'  Monographs of the Society for Research in Child Development, 78(5), i–150.
#'  
#'  Burchinal, M., Howes, C., Pianta, R., Bryant, D., Early, D., Clifford, R., & Barbarin, O. (2008). 
#'  Predicting Child Outcomes at the End of Kindergarten from the Quality of Pre-Kindergarten Teacher–Child Interactions and 
#'  Instruction. Applied Developmental Science, 12(3), 140–153. https://doi.org/10.1080/10888690802199418
#'  
#'@keywords datasets
"sim_data_wide"


#' Wide complete data (binary exposure)
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). These data contain economic strain (ESEATA1) as a binary variable. 
#'
#' @name sim_data_wide_bin.rda
#' @docType data
#' @format A data frame
#'
"sim_data_wide_bin"


#' Wide data with missingness (continuous exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package. 
#' These data contain economic strain (ESEATA1) as a continuously distributed variable. 
#'
#' @name sim_data_wide_miss.rda
#' @docType data
#' @format A data frame
#'
"sim_data_wide_miss"


#' Wide data with missingness (binary exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package. 
#' These data contain economic strain (ESEATA1) as a binary variable. 
#'
#' @name sim_data_wide_miss_bin.rda
#' @docType data
#' @format A data frame
#'
"sim_data_wide_miss_bin"


#' Long data with missingness (continuous exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package. 
#' These data contain economic strain (ESEATA1) as a continuously distributed variable. 
#'
#' @name sim_data_long_miss.rda
#' @docType data
#' @format A data frame
#'
"sim_data_long_miss"


#' Long data with missingness (binary exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package. 
#' These data contain economic strain (ESEATA1) as a binary variable. 
#'
#' @name sim_data_long_miss_bin.rda
#' @docType data
#' @format A data frame
#'
"sim_data_long_miss_bin"


#' Wide data imputed with mice (continuous exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package before 
#' imputing with the mice package. These data contain economic strain (ESEATA1) as a continuously distributed variable. 
#'
#' @name sim_data_mice.rda
#' @docType data
#' @format A mice object
#'
"sim_data_mice"


#' Wide data imputed and read in (continuous exposure)
#'
#' These data are simulated based on data from the Family Life Project (FLP), a longitudinal study following 1,292 families 
#' representative of two geographic areas (three counties in North Carolina and three counties in Pennsylvania) with high rural 
#' child poverty (Vernon-Feagans et al., 2013; Burchinal et al., 2008). MAR missingness has been added using the missMethods package before 
#' imputing with the mice package and reading in each imputed dataset. These data contain economic strain (ESEATA1) as a continuously 
#' distributed variable. 
#'
#' @name sim_data_imp_list.rda
#' @docType data
#' @format A list of data frames
#'
"sim_data_imp_list"

