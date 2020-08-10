\name{GiG.df.GenomicFeatures}
\alias{GiG.df.GenomicFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Genomic features and G4-iM Grinder's results (PQS and PiMS) collocation function
}
\description{
GiG.df.GenomicFeatures function is suitable to determine the genomic features that share their location with (and hence may be affected by) GiG's PQS and PiMS results. It employs the online database connector package 'biomartr' to retrieve the genomic annotations file for the sequence, with which then to match positions. The function returns a data frame, each row being a PQS or PiMS-hit per feature. Several feature characteristics and ID's are also returned to the user, and linked to their quadruplex-hits through the row.id variable identifier.
}
\usage{
GiG.df.GenomicFeatures(df, org, db, NumRow = NA, Feature = NA, sep = ";")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
data.frame, G4-iM Grinder M2a or M3a result data frame (PQSM2a or PQSM3a) of the analysis of the sequence with which colocate the biological features of the genome.
}
  \item{org}{
character, organism identificator used by "biomartr" package for the retrieval of the GFF annotations files. As stated in the original package: the character can specify the scientific name: e.g. organism = "Homo sapiens", the database specific accession identifier: e.g. organism = "GCF_000001405.37" (= NCBI RefSeq identifier for Homo sapiens) or the taxonomic identifier from NCBI Taxonomy: e.g. organism = "9606" (= taxid of Homo sapiens)
}
  \item{db}{
character, database used by "biomartr" package to download the GFF file using the "biomartr" package. Currently supports "refseq" or "genbank".
}
\item{reference}{
logical, used by "biomartr" package to indicate whether or not a genome shall be downloaded if it isn't marked in the database as either a reference genome or a representative genome. The factory-fresh is FALSE.
}
  \item{path}{
character, used by "biomartr" package specifying the location (a folder) in which the corresponding annotation file shall be stored. The factory-fresh is path = file.path("_ncbi_downloads","genomes").
}
  \item{NumRow}{
integer or vector of integers, the rows with the positions of the PQS or PiMS in the genome to colocate with the GFF file. It can be left as \code{NA} for the analysis of the entire dataframe. The factory-fresh default is \code{NA}.
}
  \item{Feature}{
character or vector of character, the genomic features to colocate with G4-iM Grinder's results (e.g. "CDS" or c("gene", "UTR", "telomere")). It can be left as \code{NA} for the analysis of all genomic features. The factory-fresh default is \code{NA}.
}
  \item{sep}{
character, separator character used to split information from the GFF attribute column. The factory-fresh default is \code{";"}.
}

}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result of GiG.df.GenomicFeatures is a data.frame with the collocated hits between PQS or PiMS detected with G4-iM Grinder and the genomic features. Each row is a hit between a result and a genomic feature.
}
\section{Column Meanings}{
\emph{row.id}{: row number of PQS or PiMS in the input df which has been colocated with a genomic feature. Integer used as ID number with which link the df results with the GiG.df.GenomicFeatures results.}

\emph{Start}{: start position of the PQS or PiMS which has had a hit with a genomic feature.}

\emph{Finish}{: finish position of the PQS or PiMS which has had a hit with a genomic feature.}

\emph{Score}{: score of the PQS or PiMS which has had a hit with a genomic feature.}

\emph{Conf.Quad.Seqs}{: name and counts of known-to-form-quadruplex sequences within a PQS or PiMS which has had a hit with a genomic feature.}

\emph{Strand}{: strand position of the PQS or PiMS which has had a hit with a genomic feature.}

\emph{type}{: biological feature type which has been colocated with the PQS or PiMS.}

\emph{Attributes.X}{: attributes of the biological feature which has been colocated with the PQS or PiMS.}

}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Efres Belmonte-Reche
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

  # Downloading a sequence from refseq via biomartr package
  require(G4iMGrinder)
  require(biomartr)

  #Ebola virus genome identificator for the refseq database
  org <- "GCF_000855585.1"

  Sequence <- toString(biomartr::read_genome(biomartr::getGenome(db = "refseq", organism = "GCF_000855585.1", reference = F)))
  Name <- names(biomartr::read_genome(biomartr::getGenome(db = "refseq", organism = "GCF_000855585.1", reference = F)))

  # Creating the G4iMGrinder Results for the DNA search of G-quadruplex.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence)

  # Finding biological features that colocate with the PQS results (overlapping size-restricted search).
  Rs.Features <- GiG.df.GenomicFeatures(df = Rs$PQSM2a, org = org, db = "refseq")

  # Finding biological features that colocate with the Potential Higher Order Quadruplex Sequence (PHOQS) results (non-overlapping size-unrestricted search).
  Rs.Features <- GiG.df.GenomicFeatures(df = Rs$PQSM3a, org = org, db = "refseq")

  ### Please consider that the analysis cannot be done for the PQSM2b and PQSM3b dataframes as they do not have starting nor finishing positions.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
