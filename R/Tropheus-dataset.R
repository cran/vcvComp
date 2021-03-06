#' Tropheus dataset
#'
#' A data frame of 723 observations of 57 variables extracted from a freely available dataset,
#' downloaded from the Dryad digital repository (\url{https://doi.org/10.5061/dryad.fc02f}).
#' The observations correspond to cichlid fishes of the species \emph{Tropheus moorii}
#' (color morphs 'Kaiser' and 'Kirschfleck') and \emph{T. polli} collected from eight locations
#' of Lake Tanganyika (Kerschbaumer et al., 2014).
#' The main numerical variables provided are the 2D Cartesian coordinates of 19 landmarks
#' quantifying the external body morphology of adult fishes
#' and the genotypes for 6 microsatellite markers.
#'
#' \itemize{
#' \item \strong{List_TropheusData_ID} {Specimen ID}
#' \item \strong{Extractionnr.} {Extraction number for genomic DNA}
#' \item \strong{G} {Group number}
#' \item \strong{POP.ID} {Population Id}
#' \item \strong{Sex} {Sex}
#' \item \strong{Allo.Symp} {Allopatric or sympatric population}
#' \item \strong{X1 ... Y19} {Cartesian coordinates of 19 landmarks}
#' \item \strong{Pzep3_1 ... UME003_2} {Genotype for 6 microsatellite markers}
#' }
#'
#' @references Kerschbaumer M, Mitteroecker P, Sturmbauer C (2014)
#' Evolution of body shape in sympatric versus non-sympatric Tropheus populations of Lake Tanganyika. \emph{Heredity 112(2)}: 89–98. \url{https://doi.org/10.1038/hdy.2013.78}
#' @references Kerschbaumer M, Mitteroecker P, Sturmbauer C (2013)
#' Data from: Evolution of body shape in sympatric versus non-sympatric Tropheus populations of Lake Tanganyika. \emph{Dryad Digital Repository}. \url{https://doi.org/10.5061/dryad.fc02f}
#'
#' @name Tropheus
#' @usage data(Tropheus)
#' @format A data frame with 723 rows and 57 variables
#' @docType data
NULL
