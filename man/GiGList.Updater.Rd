\name{GiGList.Updater}
\alias{GiGList.Updater}
\title{
Analysis Updater for a G4-iM Grinder List (GiGList)
}
\description{
A function that can be used to update the results of G4-iM Grinder List with different a variable configuration.
Allows the complete analytic potential on existing results (except cGcC).
Allows detection of opposite base-pair quadruplex structures on the opposite strand, avoiding all together a new quadruplex search [it is more efficient than a new G4-iM Grinder analysis].
}
\usage{
G4iMGrinder::GiGList.Updater(GiG.List, RunComposition = "G", LoopSeq = c("G", "T", "A", "C"), G4hunter = TRUE, PQSfinder = TRUE, Bt = 14, Pb = 17, Fm = 3, Em = 1, Ts = 4, Et = 1, Is = -19, Ei = 1, Ls = -16, ET = 1, WeightParameters = c(50, 50, 0), FreqWeight = 0, KnownQuadruplex = TRUE, KnownNOTQuadruplex = FALSE, NCores = 1 )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GiGList}{
list, G4-iM Grinder results in the form of a GiGList
}
  \item{ChangeRunComposition}{
logical, allows the detection of opposite base-pair quadruplex structures. Predefined values is FALSE.
}
\item{...}{
Rest of items can be found in the G4iMGrinder documentation. These values are all set to predefined values as in G4-iM Grinder
}
}
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
The result of the function is a GiGList with the new analysis.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

  # Creating the G4iMGrinder Results for the DNA search of G-quadruplex.
  Rs.G4 <- G4iMGrinder(Name = Name, Sequence = Sequence)

  # Finding the potential i-Motifs using the previous analysis,  and including known NOT quadruplex structures (by default OFF).
  Rs.iM <- GiGList.Updater(Rs.G4, ChangeRunComposition = TRUE, KnownNotQuadruplex = TRUE)

  # Changing score weights and frequency importance in calculations of Final Score
  Rs.iMScores <- GiGList.Updater(Rs.iM, WeightParameters = c(100, 50, 0), FreqWeight = 50)
}
