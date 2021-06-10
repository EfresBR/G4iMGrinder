\name{GiG.df.GenomicFeatures}
\alias{GiG.df.GenomicFeatures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Co-localization function of G4-iM Grinder's results with the annotation features
}
\description{
The function is suitable for determining the genomic features that share their location with (and hence may be affected by) GiG's PQS and PiMS results.
 It requires of an annotation file for the sequence, with which then will match positions.
 The function returns a data frame (within a tibble) of all the matches found for the input sequences with all the information of the genomic feature hit. }
\usage{
GiG.df.GenomicFeatures(GiG.df, , db, NumRow = NA, Feature = NA, sep = ";")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GiG.df}{
data.frame, G4-iM Grinder M2a or M3a result data frame (PQSM2a or PQSM3a). Usually within a GiG.List.
}
  \item{GFF}{
data.frame, annotation file as a data.frame. Must have "Start", "Finish", "type" and "attribute" columns.
}
  \item{NumRow}{
integer or vector of integers, the rows of the GiG.df to colocate with the GFF file. It can be left as \code{NA} for the analysis of the entire dataframe. The factory-fresh default is \code{NA}.
}
  \item{Feature}{
character or vector of characters, the genomic features to colocate with G4-iM Grinder's results (e.g. "CDS" or c("gene", "UTR", "telomere")). It can be left as \code{NA} for the colocalization on all genomic features. The factory-fresh default is \code{NA}.
}
  \item{sep}{
character, separator character used to split information from the GFF attribute column. The factory-fresh default is \code{";"} (for future upgrades).
}

}
\details{
If both GiG.df and GFF files have a strand column ("Strand" or "strand" can be used), the results will be matched also by strand. Else, strands will not be used to match the results and a warning will be displayed regarding the circumstance.
Please make sure the annotation file supplied matches the sequence used in the analysis.
G4-iM Grinder reverse complements the sequence to generate the supplementary strand. If the annotation files is only the complement, the positions of one of the files has first to be adapated (see below).


%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function creates two new columns for the GiG.df.

nFeatures is the total number of biological features co-located with the results.

Features is a dataframe with the actual biological features that co-locate with the results.


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

  # Getting the annotatation file of the genome
  GFF <- biomartr::getGFF(organism = "GCF_000756205.1", reference = F)
  GFF <- biomartr::read_gff(GFF)

  # Finding biological features that colocate with the PQS results (overlapping size-restricted search).
  Rs.Features <- GiG.df.GenomicFeatures(df = Rs$PQSM2a, GFF= GFF )

  # Finding biological features that colocate with the Potential Higher Order Quadruplex Sequence (PHOQS) results (non-overlapping size-unrestricted search).
  Rs.Features <- GiG.df.GenomicFeatures(df = Rs$PQSM3a, org = org, db = "refseq")



  # To change the direction (regarding the PQS positions) of the supplementary strand in G4-iM Grinder's results.
  # Get results dataframe (Here for example, Method 2 (M2)- overlapping and size-restricted search). 
  M2 <- Rs$PQSM2a 
  #Length of Sequence
  LSeq <- nchar(Sequence)
  #Dataframe of old positions.
  Supdf <- data.frame(Start = M2$Start[M2$Strand =="-"], Finish = M2$Finish[M2$Strand =="-"])
  # invert the positions (from reverse complement, to complement).
  M2$Start[M2$Strand =="-"] <- LSeq - Supdf$Finish 
  M2$Finish[M2$Strand =="-"] <- LSeq - Supdf$Start 
  # Position are inverted. Complementary is now in anti-sense direction (complement). 
  # Then apply GiG.df.GenomicFeatures function.

  }
	
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
