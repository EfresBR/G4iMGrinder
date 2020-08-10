\name{GiG.DB.rda}
\alias{GiG.DB.rda}
\docType{data}
\title{
Quadruplex database of G4-iM Grinder.
}
\description{
Database of DNA and RNA, G4s and i-Motifs, forming and non-forming structures and their characteristics. Database is divided in four data frames.

1. GiG.DB.BioInformatics, contains the sequences to be found by G4-iM Grinder and information.

2. GiG.DB.Refs, contains the sequences' bibliographical documentation.

3. GiG.DB.BioPhysical, contains the sequences' biophysical-related information.

4. Version, contains the version information of the database.

}
\usage{
data(GiG.DB$GiG.DB.BioInformatics)
data(GiG.DB$GiG.DB.Refs)
data(GiG.DB$GiG.DB.BioPhysical)
data(GiG.DB$Version)
}
\format{...}


\section{1. GiG.DB$GiG.DB.BioInformatics}{

  \describe{
    \item{\code{Quadruplex}}{logical, does the sequence form a quadruplex. TRUE for yes, FALSE for no.}
    \item{\code{Genome}}{factor, genome type (DNA or RNA)}
    \item{\code{Nucleotide}}{factor, run-nucleotide composition (G or C)}
    \item{\code{Name}}{character, the name of the quadruplex/sequence. Must be unique in Bioinformatics dataframe. Must have entry in Refs data frame.}
    \item{\code{Sequence}}{character, the nucleotide arrangement of the sequence. Must be unique in Bioinformatics dataframe.}
    \item{\code{Length}}{integer, the length in nucleotides of the quadruplex sequence}
    \item{\code{Tm}}{integer, Number of biophysical entries within the GiG.DB$GiG.DB.BioPhysical dataframe in relationship with the sequence.}
}}


\section{2. GiG.DB$GiG.DB.Refs}{

  \describe{

    \item{\code{Name}}{character, the name of the quadruplex/sequence. Used as ID to link with Bioinformatics data frame. Must be unique in Refs data frame. Must have an equal entry in Bioinformatics dataframe.}
    \item{\code{DOI}}{character, DOI code of article where sequence was published.}
    \item{\code{Pubmed}}{character, PubMedID code of article where sequence was published.}
    \item{\code{comments}}{character, Other information and comments of article where sequence was published.}
}}


\section{3. GiG.DB$GiG.DB.BioPhysical}{

       \describe{
          \item{\code{Name}}{character, the name of the quadruplex/sequence. Used as ID to link with Bioinformatics data frame. Must have an equal entry in Bioinformatics and Refs dataframe.}
          \item{\code{Tm}}{integer, experimental melting temperature (Tm) of the entry in degrees.}
          \item{\code{pH}}{integer, experimental pH of the entry.}
          \item{\code{CC.Sequence.uM}}{integer, experimental micromolar concentration of the DNA or RNA sequence in solution of the entry.}
          \item{\code{CC.K.mM}}{integer, experimental milimolar concentration of potassium ions in solution of the entry.}
          \item{\code{CC.Na.mM}}{integer, experimental milimolar concentration of sodium ions in solution of the entry.}
          \item{\code{Topology}}{factor, found experimental topology of the entry.}
}}

\section{4. GiG.DB$Version}{
\describe{
          \item{\code{Name}}{factor, data frames within GiG.DB.}
          \item{\code{Tm}}{integer, version of data frames.}
          }}

\details{

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
\keyword{datasets,G4,DNA, RNA, i-Motifs, melting temperature}
