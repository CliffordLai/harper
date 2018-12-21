#' Microarray time series experiment for yeast cell cycle from alpha experiment
#'
#' 6,178 yeast genes expression measures (log-ratios) with 
#' series length 18 from the alpha factor experiment.
#'
#' @usage data(alpha)
#' 
#' @format Matrix with 6178 rows and 18 columns.
#' Some missing data.
#' Rows and columns are labelled.
#' - attr(*, "dimnames")=List of 2
#' ..$ : chr [1:6178] "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
#' ..$ : chr [1:18] "alpha0" "alpha7" "alpha14" "alpha21" ...
#' 
#' @references  Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R.,
#'  Anders, K., Eisen, M. B., ... & Futcher, B. (1998). 
#'  Comprehensive identification of cell cycle-regulated genes of 
#'  the yeast Saccharomyces cerevisiae by microarray hybridization. 
#'  Molecular biology of the cell, 9(12), 3273-3297.
#' 
#' Dudoit S (2016). yeastCC: Spellman et al. (1998) and 
#' Pramila/Breeden (2006) yeast cell cycle microarray data. 
#' R package version 1.12.0.
#' 
#' @source The data is extracted from the ExpressionSet of 
#' the R package \code{yeastCC}.
#' 
#' @examples 
#' data(alpha)
#' qqnorm(colMeans(alpha, na.rm=TRUE))
#' qqnorm(rowMeans(alpha, na.rm=TRUE))
"alpha"


#' Microarray time series experiment for yeast cell cycle from cdc15 experiment
#'
#' 6,178 yeast genes expression measures (log-ratios) with 
#' series length 24 from the cdc15 experiment.
#'
#' @usage data(cdc15)
#'
#' @format Matrix with 6178 rows and 24 columns.
#' Some missing data.
#' Rows and columns are labelled.
#' - attr(*, "dimnames")=List of 2
#' ..$ : chr [1:6178] "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
#' ..$ : chr [1:24] "cdc15.10" "cdc15.30" "cdc15.50" "cdc15.70" ...
#' 
#' @references  Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, 
#' V. R., Anders, K., Eisen, M. B., ... & Futcher, B. (1998). 
#' Comprehensive identification of cell cycle-regulated genes of 
#' the yeast Saccharomyces cerevisiae by microarray hybridization. 
#' Molecular biology of the cell, 9(12), 3273-3297.
#' 
#' Dudoit S (2016). yeastCC: Spellman et al. (1998) and 
#' Pramila/Breeden (2006) yeast cell cycle microarray data. 
#' R package version 1.12.0.
#' 
#' @source The data is extracted from the ExpressionSet of 
#' the R package \code{yeastCC}.
#' 
#' @examples 
#' data(cdc15)
#' qqnorm(colMeans(cdc15, na.rm=TRUE))
#' qqnorm(rowMeans(cdc15, na.rm=TRUE))
"cdc15"


#' Microarray time series experiment for yeast cell cycle from cdc28 experiment
#'
#' 6,178 yeast genes expression measures (log-ratios) with series 
#' length 17 from the cdc28 experiment.
#'
#' @usage data(cdc28)
#'
#' @format Matrix with 6178 rows and 17 columns.
#' Some missing data.
#' Rows and columns are labelled.
#' - attr(*, "dimnames")=List of 2
#' ..$ : chr [1:6178] "YAL001C" "YAL002W" "YAL003W" "YAL004W" ...
#' ..$ : chr [1:17] "cdc28.0" "cdc28.10" "cdc28.20" "cdc28.30" ...
#' 
#' @references  Spellman, P. T., Sherlock, G., Zhang, M. Q., Iyer, V. R.,
#'  Anders, K., Eisen, M. B., ... & Futcher, B. (1998). 
#'  Comprehensive identification of cell cycle-regulated genes of 
#'  the yeast Saccharomyces cerevisiae by microarray hybridization. 
#'  Molecular biology of the cell, 9(12), 3273-3297.
#' 
#' Dudoit S (2016). yeastCC: Spellman et al. (1998) and 
#' Pramila/Breeden (2006) yeast cell cycle microarray data. 
#' R package version 1.12.0.
#' 
#' @source The data is extracted from the ExpressionSet of 
#' the R package \code{yeastCC}.
#' 
#' @examples 
#' data(cdc28)
#' qqnorm(colMeans(cdc28, na.rm=TRUE))
#' qqnorm(rowMeans(cdc28, na.rm=TRUE))
"cdc28"


#' Microarray time series experiment for Caulobacter crescentus bacterial cell cycle
#'
#' In this microarray experiment there are 3062 genes measured every 1 hour.
#' There are 19 is missing gene labels and these have been given 
#' labels ORFna1,...,ORFna19.
#' There 310 with duplicate labels.
#' Of these duplicate labels, 295 are duplicated twice, 12 are duplicated 3 times
#' and 3 are duplicated 4 times.
#' Duplicate labels are renamed ORF... to ORF...a and ORF...b etc.
#'
#' @usage data(Cc)
#'
#' @format Matrix with 3062 rows and 11 columns.
#' Some missing data.
#' Rows and columns are labelled.
#' - attr(*, "dimnames")=List of 2
#' ..$ : chr [1:3062] "ORF06244a" "ORF03152a" "ORF03156a" "ORF03161a" ...
#' ..$ : chr [1:11] "1" "2" "3" "4" ...
#' 
#' @details Gene expression from synchronized cultures of the bacterium 
#' Caulobacter crescentus (Laub et al., 2000).
#' (Laub et al., 2000) identified 553 genes whose
#' messenger RNA levels varied as a function of the cell cycle but their
#' statistical analysis was not very sophisticated and they probably 
#' identified too many genes.
#' Wichert et al. (2004) found that 44 genes were found which displayed
#' periodicity based on the Fisher's g-test using a FDR with q=0.05. 
#' 
#' @references  Laub, M.T., McAdams,H.H., Feldblyum,T., 
#' Fraser,C.M. and Shapiro,L. (2000) Global analysis of the genetic 
#' network controlling a bacterial cell cycle Science,
#' 290, 2144-2148.
#' 
#' Wichert,S., Fokianos K. and Strimmer K. (2004) 
#' Identifying periodically expressed
#' transcrips in microarray time series data. Bioinformatics, 18, 5-20.
#' 
#' @examples 
#' data(Cc)
#' qqnorm(colMeans(Cc, na.rm=TRUE))
#' qqnorm(rowMeans(Cc, na.rm=TRUE))
"Cc"


#' Benchmark set B1
#'
#' List for yeast genes which are \bold{most} likely to be periodic 
#' (the benchmark set 1 in de Lichtenberg et al. (2005)).
#'
#' @usage data(B1)
#'
#' @format A vector containg 113 genes' names.
#' 
#' @details A total of 113 genes previously identified as periodically
#' expressed in small-scale experiments. The set encompasses the
#' 104 genes used by Spellman et al. (1998) and nine genes added
#' by Johansson et al. (2003).
#' 
#' @references De Lichtenberg, U., Jensen, L. J., Fausboll, A., 
#' Jensen, T. S., Bork, P.,& Brunak, S. (2005). 
#' Comparison of computational methods for the identification 
#' of cell cycle-regulated genes. Bioinformatics, 21(7), 1164-1171.
#' 
#' @source The raw data can be downloaded from 
#' \url{http://www.cbs.dtu.dk/cellcycle/yeast_benchmark/benchmark.php}.
#' 
#' @examples 
#' data(alpha)
#' data(B1)
#' alphaB1 <- alpha[rownames(alpha) \%in\% B1, ]
#' 
"B1"


#' Benchmark set B2
#'
#' List for yeast genes which are \bold{most} likely to be periodic 
#' (the benchmark set 2 in de Lichtenberg et al. (2005)).
#'
#' @usage data(B2)
#' 
#' @format A vector containg 352 genes' names.
#' 
#' @details Genes whose promoters were bound (P-value below 0.01) by at
#' least one of nine known cell cycle transcription factors in both
#' of the Chromatin IP studies by Simon et al. (2001) and Lee et al.
#' (2002). To obtain a benchmark set that is independent of B1, we
#' removed all genes included in B1 (50). The resulting benchmark
#' set, B2, consists of 352 genes of which many should be expected
#' to be cell cycle regulated, since their promoters are associated
#' with known stage-specific cell cycle transcription factors.
#' 
#' @references De Lichtenberg, U., Jensen, L. J., Fausboll, A., 
#' Jensen, T. S., Bork, P.,& Brunak, S. (2005). 
#' Comparison of computational methods for the identification 
#' of cell cycle-regulated genes. Bioinformatics, 21(7), 1164-1171.
#' 
#' @source The raw data can be downloaded from 
#' \url{http://www.cbs.dtu.dk/cellcycle/yeast_benchmark/benchmark.php}.
#' 
#' @examples 
#' data(alpha)
#' data(B2)
#' alphaB2 <- alpha[rownames(alpha)  \%in\% B2,]
#' 
"B2"


#' Benchmark set B3
#'
#' List for yeast genes which are \bold{less} likely to be periodic
#'(the benchmark set 3 in de Lichtenberg et al. (2005)).
#'
#' @usage data(B3)
#'
#' @format A vector containg 518 genes' names.
#' 
#' @details Genes annotated in MIPS (Mewes et al., 2002) as 'cell cycle
#' and DNA processing'. From these, we removed genes annotated
#' specifically as 'meiosis' and genes included in B1 (67), leaving
#' 518 genes. As a large number of genes involved in the cell cycle
#' are not subject to transcriptional regulation (not periodic), and
#' because B1 was explicitly removed, a relatively small fraction
#' of these genes should be expected to be periodically expressed.
#' 
#' @references De Lichtenberg, U., Jensen, L. J., Fausboll, A., 
#' Jensen, T. S., Bork, P.,& Brunak, S. (2005). 
#' Comparison of computational methods for the identification 
#' of cell cycle-regulated genes. Bioinformatics, 21(7), 1164-1171.
#' 
#' @source The raw data can be downloaded from 
#' \url{http://www.cbs.dtu.dk/cellcycle/yeast_benchmark/benchmark.php}.
#' 
#' @examples 
#' data(alpha)
#' data(B3)
#' alphaB3 <- alpha[rownames(alpha)  \%in\% B3,]
#' 
"B3"


#' The result for the RSR method with respect to 
#' the normal likelihood ratio test
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the normal likelihood ratio test.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^6 with 1000 replications and the
#' interpolated series sizes are c(8:19,seq(20,50,2),seq(55,100,5)).
#'
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tableRegLs"


#' The result for the RSR method with respect to 
#' the laplace likelihood ratio test (Even series length)
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the laplace likelihood ratio test.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^5 with 100 replications and the
#' interpolated series sizes are c(seq(8,50,2),seq(54,98,4),100).
#' 
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tableRegL1LaplaceEven"

#' The result for the RSR method with respect to 
#' the laplace likelihood ratio test (Odd series length)
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the laplace likelihood ratio test.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^5 with 100 replications and the
#' interpolated series sizes are c(seq(9,51,2),seq(55,99,4)).
#' 
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tableRegL1LaplaceOdd"


#' The result for the RSR method with respect to 
#' the robust g test
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the robust g test.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^5 with 200 replications and the
#' interpolated series sizes are c(8:19,seq(20,50,2),seq(55,100,5)).
#' 
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions in 
#' econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tableRgAn"


#' The result for the RSR method with respect to 
#' the Fisher'g test (Even series length)
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the Fisher's g test 
#' when series length is even.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^6 with 1000 replications and the
#' interpolated series sizes are c(seq(8,50,2),seq(54,98,4),100).
#' 
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions in 
#' econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tablegEven"


#' The result for the RSR method with respect to 
#' the Fisher'g test (Odd series length)
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the Fisher's g test 
#' when series length is odd.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^6 with 1000 replications and the
#' interpolated series sizes are c(seq(9,51,2),seq(55,99,4)).
#' 
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions in 
#' econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tablegOdd"

#' The result for the RSR method with respect to 
#' the normal likelihood ratio test to the rank transformed series
#'
#' This data set is used to build up the response surface regressions 
#' in order to obtain the p-values for the normal likelihood ratio test.
#' The simulation set-up of the response surface regression method is 
#' i.i.d. statistics of length 10^6 with 1000 replications and the
#' interpolated series sizes are c(8:19,seq(20,50,2),seq(55,100,5)).
#'
#' @format A list object with the following components:
#' 
#' \code{model_matrix}: A matrix with the estimated coefficents 
#'                      for the Response Surface Regressions (one row per quantile). 
#' 
#' \code{nc}: a vectpr of the lengths of the simulated series.
#' 
#' \code{pc}: A vector containing all probabilities at which  
#'            the quantiles are estimated.
#' 
#' \code{r}: the rate in the basis function g(n)=(1/n)^r.
#' 
#' \code{q}: the number of parameters (except the intercept) used in 
#'           the response surface regression, the intercept term is removed if negative.
#' 
#' @references MacKinnon, James (2001) : 
#' Computing numerical distribution functions 
#' in econometrics, Queen's Economics Department Working Paper, No. 1037.
#' 
#' @keywords internal
#' 
"tableRegLsRank"