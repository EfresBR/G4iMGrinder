\name{G4iMGrinder}
\alias{G4iMGrinder}
\title{
Detect and analyze possible G-quadruplex, i-Motifs and their higher order structures in DNA and RNA sequences.
}
\description{
A function to detect and analyse quadruplex sequences in a genome. G4iM Grinder can be applied as a tool for possible G-quadruplex, i-Motif and higher-order structure identification, characterization and punctuation as the probability of in vitro formation and biological relevance. The search algorithm is highly configurable in all of the process steps.
}
\references{
Belmonte-Reche,E. and Morales,J.C. (2019) G4-iM Grinder: when size and frequency matter. G-Quadruplex, i-Motif and higher order structure search and analysis tool. NAR Genomics and Bioinformatics, 2.
}
\author{
Efres Belmonte-Reche
}
\note{
M1 is Method 1. M2 is Method 2. M3 is method 3.
}
\usage{
G4iMGrinder(Name, Sequence, DNA = TRUE, Complementary = TRUE, RunComposition = "G", BulgeSize = 1, MaxIL = 3, MaxRunSize = 5, MinRunSize = 3, MinNRuns = 4,
MaxNRuns = 0, MaxPQSSize = 33, MinPQSSize = 15, MaxLoopSize = 10, MinLoopSize = 0, LoopSeq = c("G", "T", "A", "C"), Method2 = TRUE, Method3 = TRUE, G4hunter = TRUE,
cGcC = FALSE, PQSfinder = TRUE, Bt = 14, Pb = 17, Fm = 3, Em = 1, Ts = 4, Et = 1, Is = -19, Ei = 1, Ls = -16, ET = 1, WeightParameters = c(50, 50, 0),
FreqWeight = 0, KnownQuadruplex = TRUE, KnownNOTQuadruplex = FALSE,  RunFormula = FALSE, NCores = 1, Verborrea = TRUE)
}
\arguments{
  \item{Name}{
character, name of the DNA or RNA sequence to grind.
}
  \item{Sequence}{
character, DNA or RNA sequence to grind composed of the nucleotide arrangement.
}
  \item{DNA}{
logical, controls if the sequence is DNA or RNA. The factory-fresh default is \code{TRUE} assuming the sequence is DNA.
}
  \item{Complementary}{
logical, controls if the Complementary strand should be created and analyzed. The factory-fresh default is \code{TRUE}.
}
  \item{RunComposition}{
character, nucleotide that composes the runs. \code{RunComposition == "G"} for G-quadruplex and \code{RunComposition == "C"} for i-Motifs. Any other nucleotide (or letter) can be imputed. The factory-fresh default is \code{"G"}.
}
  \item{BulgeSize}{
integer, number of acceptable non-\code{RunComposition} nucleotides to exist within runs. The factory-fresh default is \code{1}. Used by M1.
}
  \item{MaxRunSize}{
integer, max. number of \code{RunComposition}-nucleotides that compose a run. The factory-fresh default is \code{5}. Used by M2.
}
  \item{MinRunSize}{
integer, min. number of \code{RunComposition}-nucleotides that compose a run. The factory-fresh default is \code{5}. Used by M1.
}
  \item{MaxLoopSize}{
integer, max. number of nucleotides that may exist between runs to assume relationship. The factory-fresh default is \code{10}. Used by M2 and M3.
}
  \item{MinLoopSize}{
integer, min. number of nucleotides that may exist between runs to assume relationship. The factory-fresh default is \code{0}. Used by M2 and M3.
}
  \item{MaxNRuns}{
integer, max. number of runs that compose a structure. The factory-fresh default is \code{0}. When \code{MaxNRuns < MinNRuns}, G4iM Grinder will evade using this \code{MaxNRuns} variable in the search algorithm, allowing looking for structures with more than the traditional 4 run. Used by M2.
}
  \item{MinNRuns}{
integer, min. number of runs that compose a structure. The factory-fresh default is \code{4}. Used by M2 and M3.
}
  \item{MaxPQSSize}{
integer, max. number of nucleotides that compose a structure. The factory-fresh default is \code{33}. Used by M2.
}
  \item{MinPQSSize}{
integer, min. number of nucleotides that compose a structure. The factory-fresh default is \code{15}. Used by Method2 and M3.
}
  \item{MaxIL}{
integer, total number of nucleotides to allow to exist in between all the \code{RunComposition}-run of a structure. The factory-fresh default is \code{3}. Used by M2.
}
  \item{Method2}{
logical, to apply method 2 (M2A) of analysis to the sequence search results. This will search for structures with defined size and runs.  Depends on variables: \code{MinNRuns}, \code{MinPQSSize}, \code{MaxNRuns} and \code{MaxPQSSize}. Will also give frequency counts of each structure detected (M2B). The factory-fresh default is \code{TRUE}
}
  \item{Method3}{
logical, to apply method 3 (M3A) to the sequence search results. Search for structures with unrestricted size and numbers of runs. Useful for searching higher forming structures. Depends on variables: \code{MinNRuns} and \code{MinPQSSize}. Will also give frequency counts of each structure detected (M3B). The factory-fresh default is \code{TRUE}
}
  \item{LoopSeq}{
character, vector that defines what nucleotide and/or nucleotide pattern to quantify in each structure detected. The factory-fresh default is \code{c("G", "C", "A", "T")} but multi-character patterns are accepted, like \code{GGG}.
}
  \item{WeightParameters}{
vector of three integers, where each of the integers are the weighted value of each possible scoring system: G4hunter, PQSfinder and cGcC (in that order). Depends on the scoring system to be \code{TRUE}, or else it will force here its value to 0. The ‘factory-fresh’ default is \code{c(0.5, 0.5, 0)}, were the final score of a structure will hence be the average between G4hunter and PQSfinder. cGcC scoring value should always be 0 as its algorithm score system is not in the 0-100 range.
}
  \item{G4hunter}{
logical, to apply G4hunter algorithm as a scoring mechanism of in vitro probability of formation. The factory-fresh default is \code{TRUE}.
}
  \item{cGcC}{
logical, to apply cGcC algorithm as a scoring mechanism of in vitro probability of formation. Only for RNA sequences. The factory-fresh default is \code{FALSE}.
}
  \item{PQSfinder}{
logical, to apply an adaptation of PQSfinder algorithm as a scoring mechanism of in vitro probability of formation. The factory-fresh default is \code{TRUE}.
}
  \item{Bt}{
integer, tetrad stacking bonus constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{14}.
}
  \item{Pb}{
integer, inter-Loop penalization constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{17}.
}
  \item{Fm}{
integer, loop length penalization constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{3}.
}
  \item{Em}{
integer, loop length exponential constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{1}.
}
  \item{Ts}{
integer, tetrad supplement constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{4}.
}
  \item{Et}{
integer, inter-Loop supplement constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{-19}.
}
  \item{Is}{
integer, loop supplement constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{-16}.
}
  \item{Ei}{
integer, tetrad exponential constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{1}.
}
  \item{Ls}{
integer, inter-Loop exponential constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{1}.
}
  \item{ET}{
integer, total formula exponential constant used for the PQSfinder adaptation calculations. The factory-fresh default is \code{1}.
}
  \item{KnownQuadruplex}{
logical, controls if G4iM Grinder should compare the results with a list of known sequences that have already been demonstrated to form in vitro. Only for \code{RunComposition = "G"} or \code{RunComposition = "C"}. The factory-fresh default is \code{TRUE}.
}
  \item{KnownNOTQuadruplex}{
logical, controls if G4iM Grinder should compare the results with a list of known sequences that have already been demonstrated to NOT form in vitro. Only for \code{RunComposition = "G"} or \code{RunComposition = "C"}. The factory-fresh default is \code{FALSE}.
}
  \item{FreqWeight}{
integer, an arbitrary constant to which calculate the importance of the structure frequency. Useful only for M2B and M3B, were frequency of the structures are calculated and a new score is computed considering structure frequency.The factory-fresh default is \code{0}.
}
  \item{RunFormula}{
logical, should the formula of the PQS be calculated. The factory-fresh default is \code{FALSE}.
}
  \item{NCores}{
integer, number of Cores to cede to the function for parallel computation. The factory-fresh default is \code{1}
}
  \item{Verborrea}{
logical, allow the function to update the user with its progress. The factory-fresh default is \code{TRUE}
}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The result of G4iM Grinder is a \code{List}
\item{Configuration}{A \code{data.frame} with the variables and configuration used by G4iM Grinder.}
\item{FunTime}{A \code{data.frame} with the time taken by each part of G4iM Grinder.}
\item{PQSM2a}{A \code{data.frame} with the results found using M2A (method 2a). Only if \code{Method2 == TRUE}.}
\item{PQSM2b}{A \code{data.frame} with the results found using M2B considering frequency (method 2b). Only if \code{Method2 == TRUE}.}
\item{PQSM3a}{A \code{data.frame} with the results found using M3A (method 3a).  Only if \code{Method3 == TRUE}.}
\item{PQSM3b}{A \code{data.frame} with the results found using M3B considering frequency (method 3b).  Only if \code{Method3 == TRUE}.}
}

\section{Column Meanings}{

\emph{Start}{: integer, start position of the sequence in the genome. Only for M2A and M3A.}

\emph{Finish}{: integer, end position of the sequence in the genome. Only for M2A and M3A.}

\emph{Freq}{: integer, sequence frequency of appearance in the genome. Only for M2B and M3B.}

\emph{Runs}{: integer, number of runs (G-runs for PQS, C-runs for PiQS) in the sequence.}

\emph{IL}{: integer, number of bulges in the sequence.}

\emph{mRun}{: numeric, average run size.}

\emph{Sequence}{: character, sequence nucleotide arrangement.}

\emph{Length}{: integer, size in nucleotides of the sequence.}

\emph{Strand}{: character, strand position of the sequence. "+" is the original and "-" is the complementary strand. Only if Complementary = TRUE.}

\emph{G4Hunter}{: numeric, sequence score by G4Hunter. Only if G4Hunter = TRUE.}

\emph{pqsfinder}{: numeric, sequence score by GiG's PQSfinder. Only if pqsfinder = TRUE.}

\emph{cGcC}{: numeric, sequence score by GiG's cGcC. Only if cGcC = TRUE.}

\emph{Score}{: numeric, sequence score combining all selected scores and sequence frequency.}

\emph{Conf.Quad.Seqs}{: character, name and times found (in parenthesis) of known-to-form quadruplexes in the sequence. DNA Known-to-form structures have asterisk (*) after the number of times detected. RNA known-to-form structures have a circumflex (^) after the times detected. Only if KnownQuadruplex = TRUE.}

\emph{Conf.NOT.Quad.Seqs}{: character, name and times found (in parenthesis) of known-NOT-to-form quadruplexes in the sequence. DNA Known-NOT-to-form structures have asterisk (*) after the number of times detected. RNA known-NOT-to-form structures have a circumflex (^) after the times detected.Only if KnownNOTQuadruplex = TRUE.}
}

\examples{
  library(G4iMGrinder)

  #Retrieving a Sequence
  {if(!require("seqinr")){install.packages("seqinr")}}
  Name <- "LmajorESTs"
  Sequence <-
    paste0(read.fasta(file = url("http://tritrypdb.org/common/downloads/release-36/Lmajor/fasta/TriTrypDB-36_Lmajor_ESTs.fasta"),
                      as.string = TRUE, legacy.mode = TRUE, seqonly = TRUE, strip.desc = TRUE, seqtype = "DNA" ), collapse = "")

  #For G-quadruplex search in DNA.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence)

  #For G-quadruplex search in RNA.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence, DNA = FALSE, cGcC = TRUE)

  #For i-Motifs search in DNA.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence, RunComposition == "C")

  #For flexible G-quadruplex search in DNA.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence, BulgeSize = 2, MaxLoopSize = 20, MaxIL = 10)

  #Visualization of Results
  View(Rs$PQSM2a)  # To view M2A (Method2a) results. Size dependent structures with overlapping.
  View(Rs$PQSM2b)  # To view M2B (Method2b) results. Method2 with frequency considerations.
  View(Rs$PQSM3a)  # To view M3A (Method3a) results. Size independent structures without overlapping.
  View(Rs$PQSM3b)  # To view M3B (Method3b) results. Method3 with frequency considerations.
}

\keyword{ datagen }
\keyword{ utilities }
