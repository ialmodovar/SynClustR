\name{Ecoli}
\alias{ecoli}
\docType{data}
\title{E. coli protein localization sites}
\description{
Protein Localization Sites

}
\usage{
data(Ecoli)
}
\format{

  Data format: 336 and 8 variables 
  
1. Sequence Name: Accession number for the SWISS-PROT database.
2. MCG: McGeoch's method for signal sequence recognition.
3. GVH: von Heijne's method for signal sequence recognition.
4. LIP: von Heijne's Signal Peptidase II consensus sequence score. Binary attribute.
5. CHG: Presence of charge on N-terminus of predicted lipoproteins. Binary attribute.
6. AAC: score of discriminant analysis of the amino acid content of outer membrane and periplasmic proteins.
7. ALM1: score of the ALOM membrane spanning region prediction program.
8. ALM2: score of ALOM program after excluding putative cleavable signal regions from the sequence.

Class Distribution. The class is the localization site. Please see Nakai &  Kanehisa referenced above for more details.

  cp  (cytoplasm)                                    143
  im  (inner membrane without signal sequence)        77               
  pp  (perisplasm)                                    52
  imU (inner membrane, uncleavable signal sequence)   35
  om  (outer membrane)                                20
  omL (outer membrane lipoprotein)                     5
  imL (inner membrane lipoprotein)                     2
  imS (inner membrane, cleavable signal sequence)      2


}

\references{
Paul Horton & Kenta Nakai. "A Probablistic Classification System for Predicting the Cellular Localization Sites of Proteins".Intelligent Systems in Molecular Biology, 109-115. St. Louis, USA 1996. 
}
\keyword{ecoli}