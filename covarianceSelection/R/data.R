#' 102 detected risk genes based on Satterstrom 2019
#'
#' @name validated_genes
#' @docType data
#' @format A data frame with one column
#' \describe{
#'   \item{Gene}{gene name}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source Genes based on the findings of
#' "Satterstrom, F. Kyle, et al. "Large-scale exome sequencing study implicates both developmental 
#' and functional changes in the neurobiology of autism." (2019)."
#' @keywords data
NULL

#' TADA data from De Rubeis 2014.
#' 
#' This includes 33 risk genes if you look for genes with a \code{qvalue} less than 0.1.
#'
#' @name tada
#' @docType data
#' @format A data frame with four columns
#' \describe{
#'   \item{Gene}{gene name}
#'   \item{dn.LoF}{number of measured de novo loss-of-function mutations}
#'   \item{qvalue}{q-values based on the p-values}
#'   \item{pval.TADA}{p-values based on the TADA model}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source Based on
#' "De Rubeis, Silvia, et al. "Synaptic, transcriptional and chromatin genes disrupted in autism." 
#' Nature 515.7526 (2014): 209."
#' @keywords data
NULL

#' Partitioning of the regions in the brain
#'
#' @name region_subregion
#' @docType data
#' @format A data frame with two columns
#' \describe{
#'   \item{region}{region}
#'   \item{subregion}{subregion}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The partitionings were based on
#' "Willsey, A. J., Sanders, S. J., Li, M., Dong, S., Tebbenkamp, A. T., Muhle,
#' R. A., Reilly, S. K., Lin, L., Fertuzinhos, S., and Miller, J. A. (2013).
#' Coexpression networks implicate human midfetal deep cortical projection
#' neurons in the pathogenesis of autism. Cell, 155(5):997–1007."
#' @keywords data
NULL

#' Brainspan ID information
#'
#' @name brainspan_id
#' @docType data
#' @format A data frame with 57 rows (57 ID's) and 9 columns. The ethnicity codes
#' key is not available at the present time.
#' \describe{
#'   \item{Stage}{integer to denote which stage of development the post-mortem brain is in}
#'   \item{Braincode}{ID}
#'   \item{Age}{age in PCW (post-conception weeks) or in M (months) or Y (years) after birth}
#'   \item{Days}{age but converted into days}
#'   \item{Gender}{gender}
#'   \item{Hemisphere}{L or R to denote if only the left or right hemisphere of the post-mortem
#'   brain was collected. Otherwise, both hemispheres available}
#'   \item{pH}{pH of the tissue}
#'   \item{Ethn.}{code to denote ethnicity}
#'   \item{PMI}{post-mortem interval}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The dataset is based on
#' "Kang, H. J., Kawasawa, Y. I., Cheng, F., Zhu, Y., Xu, X., Li, M., Sousa,
#' A. M., Pletikos, M., Meyer, K. A., and Sedmak, G. (2011). Spatio-temporal
#' transcriptome of the human brain. Nature, 478(7370):483–489." The raw data
#' (prior to our lab's preprocessing) is located at \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25219}.
#' @keywords data
NULL

#' Genes expressed in brain
#'
#' @name brain_expression
#' @docType data
#' @format A data frame with 22345 rows (22345 genes) and 2 columns
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{Brain_expressed}{Status of whether or not the gene is known to be expressed in the brain}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The dataset came from Bernie Devlin.
#' @keywords data 
NULL