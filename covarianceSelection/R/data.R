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

#' Gene synonyms
#'
#' @name synonyms
#' @docType data
#' @format A list
#' \describe{
#'   \item{hash}{A \code{hash} object where each key is a symbol and each value is an integer.
#'   The integer references which element in \code{syn.list} corresponds to the key.}
#'   \item{syn.list}{A list where each element (named) is a vector of synonym symbols (characters)}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The dataset is based on outputs in \url{http://biomart.genenames.org/martform/#!/default/HGNC?datasets=hgnc_gene_mart}.
#' @keywords data
NULL

#' TADA dataset
#'
#' @name tada
#' @docType data
#' @format A data frame with 18735 rows (18735 genes) and 17 columns.
#' \describe{
#'   \item{Gene}{Gene name}
#'   \item{mut.rate}{mutation rate}
#'   \item{dn.LoF}{number of denovo loss of function events}
#'   \item{case.LoF}{number of loss of function events in case}
#'   \item{ctrl.LoF}{number of loss of function events in control}
#'   \item{trans.LoF}{number of transmitted loss of function events}
#'   \item{ntrans.LoF}{number of nontransmitted loss of function events}
#'   \item{dn.mis3}{number of denovo missense 3 mutation events}
#'   \item{case.mis3}{number of missense 3 mutation events in case}
#'   \item{ctrl.mis3}{number of missense 3 mutation events in control}
#'   \item{trans.mis3}{number of transmitted missense 3 events}
#'   \item{ntrans.mis3}{number of nontransmitted missense 3 events}
#'   \item{BF.dn}{Bayes factor for denovo}
#'   \item{BF}{Bayes factor}
#'   \item{qvalue.dn}{q value for denovo}
#'   \item{qvalue}{qvalue}
#'   \item{pval.TADA}{pvalue}
#' }
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The methodology to calculate these values were based on
#' "He, X., Sanders, S. J., Liu, L., De Rubeis, S., Lim, E. T., Sutcliffe, J. S.,
#' ... & Devlin, B. (2013). Integrated model of de novo and inherited genetic
#' variants yields greater power to identify risk genes. PLoS genetics, 9(8), e1003671.",
#' but the dataset itself comes from 
#' "De Rubeis, Silvia, et al. "Synaptic, transcriptional and chromatin genes disrupted in autism." Nature 515.7526 (2014): 209."
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

#' Iossifov TADA genes
#'
#' @name iossifov
#' @docType data
#' @format A vector of characters, length 251
#' @author Kevin Lin \email{kevinl1@andrew.cmu.edu}
#' @source The partitionings were based on the data used in
#' "Iossifov, I., ORoak, B. J., Sanders, S. J., Ronemus, M., Krumm, N., Levy,
#' D., Stessman, H. A., Witherspoon, K. T., Vives, L., and Patterson, K. E. (2014).
#' The contribution of de novo coding mutations to autism spectrum disorder.
#' Nature, 515(7526):216–221.

#' @keywords data
NULL
