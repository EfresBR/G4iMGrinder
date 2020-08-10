\name{GiGList.Analysis}
\alias{GiGList.Analysis}
\title{
Summary generator for a G4-iM Grinder resulting list (GiGList)
}
\description{
To summarize the results of a GiGList into a single row. Allows easier comparison between genomes.
}
\usage{
GiGList.Analysis(GiGList, iden, ScoreMin = c(20, 40), FreqMin = 10, LengthMin = 50, Density = 100000, byDensity = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GiGList}{
List, G4iMGrinder results in the form of a GiGList to analyse.
}
  \item{iden}{
character, identification tag of the DNA or RNA origin.
}
  \item{ScoreMin}{
integer or vector of integers, score filter to apply to the results. The function will create a new column for each integer in the vector with the filtrate results that score at least the value. Values are interpreted as absolute; if the search involves i-Motif detection, it will multiply by -1 to adapt to its scoring scale.
The factory-fresh default is \code{c(20, 40)}, returning the values of sequences that have at least high probability of formation (\code{Scoremin >= 40} for PQS; \code{Scoremin <= -40} for PiMS), and at least medium probability of formation (\code{Scoremin >= 20} for PQS; \code{Scoremin <= -20} for PiMS).
}
  \item{FreqMin}{
  integer or vector of integers, frequency filter to apply to the results. The function will create a new column for each integer in the vector with the filtrate results that have a frequency of appearence in the genome of at least the value/s. The factory-fresh default is \code{10}, returning the results detected with \code{Freqmin >= 10}.
}
  \item{LengthMin}{
integer or vector of integers, length filter to apply to the results. The function will create a new column for each integer in the vector with the filtrate results that have a frequency of appearence in the genome of at least the value/s. Only for Method 3 results (M3a and M3b). The factory-fresh default is \code{50}, returning the results detected with \code{Lengthmin >= 50}.
}
\item{Density}{
integer, constant to calculate density results. The factory-fresh default is \code{100000}, returning results found per 100000 nucleotides. Only pertinent if \code{byDensity = TRUE}.
}
  \item{byDensity}{
logical, should the results be returned as a density. Calculated as Result.Density = (Density*Results)/(total genomic length). If set to \code{FALSE}, it will return counts. Density allows genomic size-indepedent comparisons. The factory-fresh default is \code{TRUE}.
}}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\references{
%% ~put references to the literature/web site here ~
}
\value{
The result of GiGList.Analysis is a one row data.frame with the summary of the G4-iM Grinder List.
}

\section{Column Meanings}{

\emph{Name}{: name of summary, given by the name of the GiGList input.}

\emph{iden}{: identification of the genome, given by iden input.}

\emph{Length}{: Length of the genomic sequence. }

\emph{SeqG}{: percentage of the genomic sequence which is G.}

\emph{SeqC}{: percentage of the genomic sequence which is C.}

\emph{nM2a}{ (column name prefix): results for Method 2a results if method is active. Returned as counts or density.}

\emph{nM2b}{ (column name prefix): results for Method 2b results if method is active. Returned as counts or density.}

\emph{nM3a}{ (column name prefix): results for Method 3a results if method is active. Returned as counts or density.}

\emph{nM3b}{ (column name prefix): results for Method 3b results if method is active. Returned as counts or density.}

\emph{.S|X|}{ (column name suffix): results for method filtered by score criteria. Within bars, the criteria. Returned as counts or density.}

\emph{.F|X|}{ (column name suffix): results for method filtered by frequency criteria. Within bars, the criteria. Returned as counts or density.}

\emph{.L|X|}{ (column name suffix): results for method filtered by length criteria. Within bars, the criteria. Returned as counts or density.}

\emph{.KTFQ}{ (column name suffix): results for method which have known-to-form-quadruplex within. Returned as counts or density.}

\emph{.KNTFQ}{ (column name suffix): results for method which have known-NOT-to-form-quadruplex within. Returned as counts or density.}

\emph{.UniqPercent}{ (column name suffix): percentage of result uniqueness per method. Calculated as: 100*[Number results with frequency of 1]/[total number of results]}

\emph{Config}{: column with analysis information of summary.}
}

\author{
Efres Belmonte-Reche
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\section{Warning }{
Column names depend on summary conditions and hence the same analysis conditions must be used in order to bind GiG.List analysis together. The analysis with different G4-iM Grinder conditions may also give rise to problems when binding results.
}
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  # Creating the G4iMGrinder Results for the DNA search of G-quadruplex.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence)

  # Creating the summary data.table
  Summary <- GiGList.Analysis(GiGList = Rs, iden = "Parasite")


}
