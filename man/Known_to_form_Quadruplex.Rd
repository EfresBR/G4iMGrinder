\name{Known_to_form_Quadruplex.rda}
\alias{Known_to_form_Quadruplex.rda}
\docType{data}
\title{
G4 forming Sequences in DNA
}
\description{
A data frame with the characteristics and refs. of in vitro-forming quadruplex sequences.
Used by G4-iM Grinder to identify results which are already known to form.
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
    \item{\code{Tm}}{integer, melting temperature as recorded in refs. pH independent for G4s}
     \item{\code{Tm.pH5}}{integer, melting temperature as recorded in refs. at pH 5, for i-Motifs}
    \item{\code{Tm.pH6}}{integer, melting temperature as recorded in refs. at pH 6, for i-Motifs}
    \item{\code{Tm.pH6.5}}{integer, melting temperature as recorded in refs. at pH 6.5, for i-Motifs}
    \item{\code{Tm.pH7}}{integer, melting temperature as recorded in refs. at pH 7, for i-Motifs}
    \item{\code{Ref}}{character, literature ref. of quadruplex}
    \item{\code{DOI}}{character, DOI number of literature ref. of quadruplex}
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
data(Known_to_form_Quadruplex)
View(Known_to_form_Quadruplex)
}
\keyword{datasets,G4,DNA}
