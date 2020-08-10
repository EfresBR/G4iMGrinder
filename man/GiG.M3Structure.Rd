\name{GiG.M3Structure}
\alias{GiG.M3Structure}
\title{
Potential Higher Order Quadruplex Sequence (PHOQS) Arrangement Analyzer for Method 3A (M3A) candidates.
}
\description{
  The purpose of this function is to locate the most probable conformation of candidates to form the PHOQS. It uses three different methods.

  HRA Conformation will select the candidate by assigning seats sequentially to the highest scoring sub units with known-to-form sequences within.

  RAH conformation will select the best conformation of candidates by locating within all possible conformations the highest mean scoring structure.

  RAnH conformation will select the best conformation of candidates by locating within all possible conformations the highest normalized mean scoring structure. The mean score of each conformation is multiplied by the percentage of the PHOQS occupied by PQS sub units.
}
\usage{
GiG.M3Structure(GiGList, M3ACandidate, MAXite)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GiGList}{
List, G4iMGrinder results in the form of a GiGList.
}
  \item{M3ACandidate}{
integer, PHOQS candidate to analyze, given as the row number it occupies in the GiGList$PQSM3a data frame.
}
  \item{MAXite}{
integer, number of iterations to repeat random seat allocation to explore all possible PHOQS arrangements. Longer and with more sub units PHOQS will require more iterations. For HoEBR1, (length: 118, number of potential PQS: 32, max.Seats: 4), 10000 iterations were enough to find all 307 arrangements.
Random allocation of seats does not contemplate vacant seats.
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
  The result of GiG.M3Structure is a list with the analysis of the selected PHOQS.

  \item{M2}{data.frame, stores the found potential sub units (located with M2A) in the location of the desired PHOQS (M3A). An error will be given if no sub units were found.}

  \item{M3}{data.frame, stores the desired PHOQS to analyze. Only 1 row is accepted.}

  \item{Potential.Arrangements}{data.frame, stores all the sub unit potential conformations that can give rise to the PHOQS. These depend on number of iterations of random seat allocation given by \code{MAXite}. Please make sure all conformations are found - See \code{MAXite}.
  }

  \item{Best.Arrangements}{list, stores the best PHOQS conformations for the three analysis, HSA conformation, RAH conformation and RAnH conformation.
}
}


\section{Column Meanings}{

\emph{RES}{: reference, of sub units (by row name) in the original GiGList$PQSM2A data frame.}

\emph{nPQS}{: integer, number of sub units that forms that configuration.}

\emph{MeanScore}{: numeric, mean score of the sub units that forms that configuration.}

\emph{PQSLenghPercent}{: numeric, percentage of the PHOQS that is involved as the sub units of that configuration.}

\emph{idenPQS}{: reference, of sub units (by row name) that form that configuration in the M2 data frame.}

\emph{nMS}{: numeric, normalized mean score of the sub units that forms that configuration taking into account the length of the PHOQS that is involved as the sub units.}
}

\examples{

# Creating the G4iMGrinder Results for the DNA search of G-quadruplex.
  Rs <- G4iMGrinder(Name = Name, Sequence = Sequence)

# Analyzing the first PHOQS found with with G4-iM Grinder
  firstPHOQS2analyze <- GiG.M3Structure(Namedf = "Rs", M3Candidate = 1, MAXite = 10000)

}
