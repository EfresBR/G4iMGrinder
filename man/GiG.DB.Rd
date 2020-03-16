\name{GiG.DB.rda}
\alias{GiG.DB.rda}
\docType{data}
\title{
Quadruplex database of G4-iM Grinder. In it, DNA and RNA Structures, G4 and i-Motif, and forming and non-forming structures are divided by their bioinformatic, biophysical and bibliographical characteristics.
}
\description{
A list that contains three data frames.
GiG.DB.BioInformatics contains the sequences to be found by G4-iM Grinder.
GiG.DB.Refs contains the sequences' bibliographical documentation.
GiG.DF.BioPhysical contains the sequences' Tm related information.
}
\usage{data(Known_to_form_Quadruplex)}
\format{
  \describe{
    \item{\code{Quadruplex}}{logical, does the sequence form a quadruplex}
    \item{\code{Genome}}{factor, examining genome type (DNA or RNA)}
    \item{\code{Nucleotide}}{factor, run-nucleotide composition (G or C)}
    \item{\code{Name}}{character, the name of the quadruplex}
    \item{\code{Sequence}}{character, the nucleotide arrangement of the sequence}
    \item{\code{Length}}{integer, the length in nucleotides of the quadruplex sequence}
    \item{\code{Tm}}{integer, melting temperature}
    \item{\code{pH}}{integer, pH value for Tm calculation}
    \item{\code{DOI}}{character, DOI identification of the ref.}
    \item{\code{Pubmed}}{character, PMID identification of the ref.}
    \item{\code{comments}}{character, other comments regarding ref.}
  }
}
\details{
\code{Name} and \code{Sequence} should be unique for each new row.

}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(GiG.DB)
View(GiG.DB)
}
\keyword{datasets,G4,DNA}
