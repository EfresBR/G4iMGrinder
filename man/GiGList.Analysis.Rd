\name{GiGList.Analysis}
\alias{GiGList.Analysis}
\title{
Summary generator for a G4-iM Grinder resulting list (GiGList)
}
\description{
To summarize the results of a GiGList into a single row. Allows easier comparison between genomes.
}
\usage{
GiGList.Analysis(GiGList, iden, ScoreMin = 40, FreqMin = 10, LengthMin = 50, Density = 1e+05, Verborrea = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GiGList}{
List, G4iMGrinder results in the form of a GiGList
}
  \item{iden}{
character, identification tag of the DNA or RNA origin.
}
  \item{Density}{
integer, constant to which calculate each method's structure density in the sequence. The factory-fresh default is \code{100000}, which hence will give the structures found per 100000 nucleotides.
}
  \item{ScoreMin}{
integer, min. score of the detected structures acceptable in a second subset of the results. The factory-fresh default is \code{40}, which means  will give the densities of the structures detected with \code{Scoremin > 40}.
}
  \item{FreqMin}{
integer, min. frequency of the detected structures acceptable in a second subset of the results. The factory-fresh default is \code{10}, which means  will give the densities of the structures detected with \code{Freqmin > 10}.
}
  \item{LengthMin}{
integer, min. length of the detected structures acceptable in a second subset of the results. The factory-fresh default is \code{50}, which means  will give the densities of the structures detected with \code{Lengthmin > 50}.
}
  \item{Verborrea}{
logical, as to let the function explain the resulting table and the meaning of each column. The factory-fresh default is \code{TRUE}.
}}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Efres Belmonte-Reche
}
\value{
The result of GiGList.Analysis is a data.frame with the G4iM Grinder List summary.

\code{Name, iden, Length, SeqG, SeqC} columns show the info of the sequence, where Length is the total length in nucleotides, and SeqG, SeqC, is the composition (as a percentage) of the sequence which is G and C respectively.

If both \code{Method2 == TRUE} and \code{Method3 == TRUE} in the G4iMGrinder results, total number of columns of the resulting df will be 27.

The \code{(nM2a)} columns show the number of results found with Method2a (M2A).

The \code{(nM2b)} columns show the number of results found with Method2b (M2B).

The \code{(nM3a)} columns show the number of results found with Method3a (M3A).

The \code{(nM3b)} columns show the number of results found with Method3b (M3B).

Suffixed (.D) columns show the density of the data. Hence, \code{nM2a.D} will show the density of all the results found by \code{nM2a}. To calculate this example, \code{nM2a.D = (nM2a * Density) / Length}.

Suffixed (.S) columns show the number of results found with the method that have a \code{FinalScore > ScoreMin}.

Suffixed (.F) columns show the number of results found with methods B(M2B or M3B) that have a \code{Frequency > FrequencyMin}.

Suffixed (.L) columns show the number of results found for method3 (M3A or M3B) that have a \code{Length > LengthMin}.

\code{(Config)} column show the configuration used to summarize the results.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  # Creating the G4iMGrinder Results for the DNA search of G-quadruplex.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence)

  # Creating the summary data.table
  Summary <- GiGList.Analysis(GiGList = Rs, iden = "Parasite")

   # Adding another summary report (with more stringent variables) to an existent result table
  Summary[,2] <- GiGList.Analysis(GiGList = Rs ,iden = "Parasite" , ScoreMin = 75, FreqMin = 50, LengthMin = 500, Density = 100000)

}
