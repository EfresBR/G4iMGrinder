####
####
##########################################################################
#################            G4-iM Grinder         #######################
##########################################################################
####    Efres Belmonte Reche            efresbr@gmail.com             ####
####    Revised                         2020-03-15                    ####
####    Version                         1.5.9                         ####
##########################################################################
##########################################################################
####
######################################################
####         ¡Santiago y Cierra, Espana!          ####
######################################################
####
####
####
####  Log
####  - G4-iM Grinder will save version of DDBB and GiG function in each Results within Configuration table.
####  - Changed KTFQS function, to lower RAM usage by dividing large datasets in 5000 chunks
####  - G4-iM Grinder will convert DNA sequences to RNA if DNA is FALSE, automatically (changing T for U). Also Reverse. Line 127-129.
#### V1.5.8
####  - Divided Known to Form Sequences DDBB into three DDBB. TmDF stores Tm values of each Sequence, REf stores the bibliographical aspects, and BioInformatic for G4-iM grinder use..
####  - Added 800 new Seqs to BioInformatic and Refs, as reversals of the original secs.
####  - Created a function to manage the DB.
####  - Updated G4-iM Grinder REf DT, Biophysical DT,  and README(!!!!). Updated subfunction for Quadruplex detection.
#### V1.5.7
#### - Fixed Bug in .PQSfinder function for when analysing quadruplexes with no bulges (BulgeSize = 0).
#### V1.5.6
#### - Fixed GiGList.Analyzer for when analysing empty GiGLists.
#### V1.5.5
#### - Fixed KnownQuadruplex finder motor to avoid error when n? of results is 1.
#### - Added safeguards to G4iMGrinder in M1A, M1B, M2A and M3A.
#### - Changed .M1A procedure to avoid errors when no results are found. Added error when no results are found.
#### - Added error on .M1B is no results are found.
#### V1.5.0
#### - Maintained mRun in M2 and M3 for posterior recalculations of PQSfinder.
#### - Changed names of complementary functions, and changed input to GiGLists.
#### - GiG.Updater: A function to update Scorering systems and known to form not Quadruplex structures for already existent results.
#### - Install instructions update. Package update. R 3.6 update
#### - Update DDBB of list of quadruplex structures to be able to form. Included 100 i-Motifs. G4Hunter and G4RNA DDBB have been read again & analyzed again - Including Known not to form Sequences.
#### - Changed var names of KnownG4 to KnownQuadruplex. Updated Var names of internal functions related.
#### - Added the detection on Known-NOT-to-form Quadruplex. Know-Not to form Quadruplex variable added. Is turned OFF unless specifically activated.
#### - Perfected The Known-To-form Quadruplex finder function (KnownQuadruplex).
#### V1.1.0
#### - Changed G4RNA name to cGcG
#### - Added function GiG.M3AStructure and documentation
#### - Fixed several bugs on M1 process regarding when resulting data is empty. Some still need to be fixed.
#### - Eliminated all comments and sounds for G4iM Grinder when verborrea = FALSE
#### V1.0.1
#### - Changed several names of varriables.
####
################################################################################
#########  G4-iM Grinder: Start of function.
################################################################################
####
G4iMGrinder <- function(Name, Sequence, DNA=TRUE, Complementary=TRUE, RunComposition="G", BulgeSize=1,
                        MaxRunSize=5, MinRunSize=3, MaxLoopSize=10, MinLoopSize=0,
                        MaxNRuns=0, MinNRuns=4, MaxPQSSize=33, MinPQSSize=15, MaxIL = 3,
                         Method2=TRUE, Method3=TRUE, LoopSeq=c("G", "T", "A", "C"),
                        WeightParameters=c(50,50,0), G4hunter=TRUE,  cGcC=FALSE, PQSfinder=TRUE,
                        Bt=14, Pb=17, Fm=3, Em=1, Ts=4, Et=1, Is=-19, Ei=1, Ls=-16, ET = 1,
                        KnownQuadruplex= TRUE, KnownNOTQuadruplex= FALSE, FreqWeight=0, RunFormula = FALSE, NCores = 1, Verborrea = TRUE){
  if(Method2 == TRUE|Method3 == TRUE){

    ### SUBFUNCTIONS OF G4iM GRINDER
    # Exported to Subfunction.R script
    #

    ######  Packages
    # Exported to Subfunction.R script
    G4iMGrinder:::.PackageLoading()

    ######  Limitations
    {
      stopifnot(is.character(Sequence), !is.null(Name), !is.null(Sequence))
      if(!is.logical(cGcC)){cGcC <- FALSE}
      if(!is.logical(KnownQuadruplex)){KnownQuadruplex<- FALSE}
      if(!is.logical(KnownNOTQuadruplex)){KnownNOTQuadruplex<- FALSE}
      if(!is.logical(G4hunter)){G4hunter <- TRUE}
      if(!is.logical(PQSfinder)){PQSfinder <- TRUE}
      if(!is.logical(Method2)){Method2 <- TRUE}
      if(!is.logical(Method3)){Method3 <- TRUE}
      if(!is.logical(DNA)){DNA <- TRUE}
      if(!is.logical(Complementary)){Complementary <- TRUE}
      if(BulgeSize > 3){ BulgeSize <- 3}
      if(BulgeSize < 0){ BulgeSize <- 0}
      if(MinRunSize < 2) {MinRunSize <- 2}
      if(MaxRunSize < MinRunSize) {MaxRunSize <- MinRunSize}
      if(MinNRuns < 2) {MinRunSize <- 2}
      if(MinLoopSize < 0) {MinPQSSize <- 0}
      if(MaxPQSSize < 10) {MaxPQSSize <- 10}
      if(MaxPQSSize < MinPQSSize) {MaxPQSSize <- MinPQSSize}
      RunComposition <- str_to_upper(RunComposition)
      LoopSeq <- str_to_upper(LoopSeq)
    if (!RunComposition %in% c("G","C")){
        KnownQuadruplex<- FALSE
        KnownNOTQuadruplex<- FALSE
        cGcC <- FALSE
        G4hunter <- FALSE
        PQSfinder <- FALSE
      }
      if (G4hunter == FALSE) {  WeightParameters[1] <- 0}
      if (PQSfinder == FALSE) {   WeightParameters[2] <- 0}
      if (cGcC == FALSE) {  WeightParameters[3] <- 0}


    }
    ######  Cores and parallel config
    {
      if (NCores > detectCores()| NCores < 1 |!is.numeric(NCores)) {
        NCores <- round(detectCores()*3/4,0)
      }
    }
    ######  Time control
    if(Verborrea == TRUE){
      cat(paste0("\nAnalysis of ", Name, " Started. Be Patient."))}
    Date <- Sys.time()
    Time <- c(0,0)
    Process <- c("initial Time","Total")
    ######  Sequence
    {
      Process[length(Process)+1] <- "Sequence Adaptations"
      Time[length(Time)+1] <- as.numeric(
        system.time({
          Sequence <- str_to_upper(Sequence)
          ifelse(test = DNA == TRUE,
             yes = Sequence <- stringr::str_replace_all(string = Sequence, pattern = "U", replacement = "T"),
             no = Sequence <- stringr::str_replace_all(string = Sequence, pattern = "T", replacement = "U"))
          if (Complementary == TRUE){
            Sequence2 <- .CompSeqFun(DNA = DNA, Sequence)}}
        )[3]
      )}
    ######    1.   G run Search
    {
      if(Verborrea == TRUE){cat("\nMethod 1A Started")}
      Process[length(Process)+1] <- "M1A"
      Time[length(Time)+1] <- as.numeric(system.time({
        PQSM1a <- .M1A(RunComposition = RunComposition, MaxRunSize= MaxRunSize, MinRunSize = MinRunSize, BulgeSize = BulgeSize, Sequence = Sequence)
        if(Complementary == TRUE){
          PQSM1aCOMP <- .M1A(RunComposition = RunComposition, MaxRunSize= MaxRunSize, MinRunSize = MinRunSize, BulgeSize = BulgeSize, Sequence = Sequence2)}
      })[3])
      if(Verborrea == TRUE){cat(
        paste0(" & Ended. ", as.numeric(nrow(PQSM1a)) + ifelse(Complementary == TRUE, as.numeric(nrow(PQSM1aCOMP)), 0), " results found"))
      }
      if(as.numeric(nrow(PQSM1a)) + ifelse(test = Complementary == TRUE, yes =  as.numeric(nrow(PQSM1aCOMP)),no =  0) == 0){
        stop("Error: M1A found no runs on sequence.")}

      if(Verborrea == TRUE){cat("\nMethod 1B Started")}
      Process[length(Process)+1] <- "M1B"
      Time[length(Time)+1] <- as.numeric(system.time({
        PQSM1b <- .M1B(MinLoopSize= MinLoopSize, MaxLoopSize = MaxLoopSize, df = PQSM1a, MinRunSize = MinRunSize, MinNRuns= MinNRuns)
        if(Complementary == TRUE){
          PQSM1bCOMP <- .M1B(MinLoopSize= MinLoopSize, MaxLoopSize = MaxLoopSize, df = PQSM1aCOMP, MinRunSize = MinRunSize, MinNRuns = MinNRuns)}
      })[3])
      if(Verborrea == TRUE){cat(
        paste0(" & Ended. ", as.numeric(nrow(PQSM1b)) + ifelse(Complementary == TRUE, as.numeric(nrow(PQSM1bCOMP)), 0), " results found"))}
      if(as.numeric(nrow(PQSM1b)) + ifelse(test = Complementary == TRUE, yes =  as.numeric(nrow(PQSM1bCOMP)), no = 0) == 0){
        stop("Error: M1B found no related runs on sequence, or they are not enough to build a quadruplex (rRuns < MinNRuns).")}
    }

    ######    2.   Method 2
    if(Method2 == TRUE){

      if(Verborrea == TRUE){cat("\nMethod 2A Started")}
      Process[length(Process)+1] <- "M2A"
      Time[length(Time)+1] <- as.numeric(system.time({
        PQSM2a <-.M2a(RunComposition = RunComposition, df= PQSM1b, NCores = NCores, MinPQSSize= MinPQSSize, MinNRuns = MinNRuns, MaxPQSSize= MaxPQSSize, MaxNRuns= MaxNRuns, PQSfinder = PQSfinder, RunFormula = RunFormula, Sequence = Sequence, MaxIL = MaxIL, Verborrea = Verborrea, BulgeSize = BulgeSize, G4hunter=G4hunter,  cGcC=cGcC)
        if(Complementary == TRUE){
          PQSM2aCOMP <-.M2a(RunComposition = RunComposition, df= PQSM1bCOMP, NCores = NCores, MinPQSSize= MinPQSSize, MinNRuns = MinNRuns, MaxPQSSize= MaxPQSSize, MaxNRuns= MaxNRuns, PQSfinder = PQSfinder, RunFormula = RunFormula, Sequence = Sequence2, MaxIL = MaxIL, Verborrea = Verborrea, BulgeSize = BulgeSize,  G4hunter=G4hunter,  cGcC=cGcC)
          if(nrow(PQSM2a)>0){
            PQSM2a$Strand <- "+"}
          if(nrow(PQSM2aCOMP)>0){
            PQSM2aCOMP$Strand <- "-"}
          if(nrow(PQSM2aCOMP)>0 & nrow(PQSM2a)>0){
            PQSM2a <- rbind(PQSM2a, PQSM2aCOMP)}
          if(nrow(PQSM2aCOMP)>0 & nrow(PQSM2a) == 0){
            PQSM2a <- PQSM2aCOMP}
          if(nrow(PQSM2aCOMP)== 0 & nrow(PQSM2a) > 0){
                PQSM2a <- PQSM2a}
          if(nrow(PQSM2aCOMP)== 0 & nrow(PQSM2a) == 0){
                PQSM2a <- data.frame()}
          rm(PQSM2aCOMP)
          }
        if(nrow(PQSM2a)>0){
          PQSM2a <- PQSM2a[order(PQSM2a$Start, decreasing = FALSE),]
          row.names(PQSM2a) <- seq(1:nrow(PQSM2a))
        }})[3])
      if(Verborrea == TRUE){cat(
        paste0(" & Ended. ", as.numeric(nrow(PQSM2a)), " results found"))
      }

      if(nrow(PQSM2a)>0){
        if(Verborrea == TRUE){cat("\nM2 - Scoring: ")}
        ##Scoring system
        #G4hunter calculus
        if (G4hunter == TRUE){
          if(Verborrea == TRUE){cat("G4Hunter (I")}
          Process[length(Process)+1] <- "M2-G4Hunter"
          Time[length(Time)+1] <- as.numeric(system.time(
            PQSM2a <- .G4Hunter(df = PQSM2a, NCores = NCores)
          )[3])
          if(Verborrea == TRUE){cat("/O).")}}
        #PQSfinder calculus
        if (PQSfinder == TRUE){
          if(Verborrea == TRUE){cat(" PQSfinder (I")}
          Process[length(Process)+1] <- "M2-PQSfinder"
          Time[length(Time)+1] <- as.numeric(system.time(
            PQSM2a <- .PQSfinderfun(PQSM2a, RunComposition, Bt=Bt, Pb=Pb, Fm=Fm, Em=Em, Ts=Ts, Et=Et, Is=Is, Ei=Ei, Ls=Ls, ET = ET, MinNRuns= MinNRuns)
          )[3])
          if(Verborrea == TRUE){cat("/O).")}}
        #cGc
        if(cGcC == TRUE){
          if(Verborrea == TRUE){cat(" cGcC (I")}
          Process[length(Process)+1] <- "M2-cGcC"
          Time[length(Time)+1] <- as.numeric(system.time(
            PQSM2a <- .cGcCfun(df = PQSM2a, RunComposition = RunComposition)
          )[3])
          if(Verborrea == TRUE){cat("/O).")}}
        #mean Score
        if((G4hunter == TRUE & PQSfinder == TRUE)|(G4hunter == TRUE & cGcC == TRUE)|(PQSfinder == TRUE & cGcC == TRUE)){
          if(Verborrea == TRUE){cat(" MeanScore (I/")}
          Process[length(Process)+1] <- "M2-meanScore"
          Time[length(Time)+1] <- as.numeric(system.time(
            PQSM2a <- .ScoreFun(df= PQSM2a, G4hunter = G4hunter, PQSfinder = PQSfinder, cGcC = cGcC, WeightParameters = WeightParameters)
          )[3])
          if(Verborrea == TRUE){cat("O).")}}

        if(Verborrea == TRUE){cat("\nQuantification Started")}
        ##Quantification
        Process[length(Process)+1] <- "M2-Quantification"
        Time[length(Time)+1] <- as.numeric(system.time(
          PQSM2a <- .LoopQuantification(df=PQSM2a, LoopSeq = LoopSeq)
        )[3])
        if(Verborrea == TRUE){cat(" & Ended. ")}

        ##Qualification
        # Known G4 Sequences
        if (((KnownQuadruplex== TRUE|KnownNOTQuadruplex == TRUE)|(KnownQuadruplex== TRUE & KnownNOTQuadruplex == TRUE)) & RunComposition == "G"|RunComposition == "C"){
          if(Verborrea == TRUE){cat("\nM2 - KnownQuadruplex Started")}
          Process[length(Process)+1] <- "M2-KnownQuadruplex"
          Time[length(Time)+1] <- as.numeric(system.time(
            PQSM2a <- .KnownQuadFun(df = PQSM2a, DNA = DNA, KnownNOTQuadruplex = KnownNOTQuadruplex, KnownQuadruplex = KnownQuadruplex, RunComposition = RunComposition)
          )[3])
          if(Verborrea == TRUE){cat(" & Ended.")}}

        if(Verborrea == TRUE){cat("\nMethod 2B Started")}
        ## Method 2B by Frequency
        Process[length(Process)+1] <- "M2-M2B"
        Time[length(Time)+1] <- as.numeric(system.time(
          PQSM2b <- .M2B(df = PQSM2a, FreqWeight = FreqWeight, RunComposition = RunComposition)
        )[3])
        if(Verborrea == TRUE){cat(" & Ended.")}} else {cat("\nNo Results found with M2.")}

    }
    ######    3.   Method 3
    if(Method3 == TRUE){

      if(Verborrea == TRUE){cat("\nMethod 3  Started")}
      Process[length(Process)+1] <- "M3"
      Time[length(Time)+1] <- as.numeric(system.time({
        PQSM3a <- .M3(df = PQSM1b, NCores = NCores, MinNRuns = MinNRuns, MinPQSSize = MinPQSSize, RunFormula = RunFormula, Sequence = Sequence, Verborrea = Verborrea, BulgeSize = BulgeSize, G4hunter=G4hunter,  cGcC=cGcC, PQSfinder = PQSfinder, RunComposition = RunComposition)
        if(Complementary == TRUE){
          PQSM3aCOMP <-.M3(df = PQSM1bCOMP, NCores = NCores, MinNRuns = MinNRuns, MinPQSSize = MinPQSSize, RunFormula = RunFormula, Sequence = Sequence2, Verborrea = Verborrea, BulgeSize = BulgeSize, G4hunter=G4hunter,  cGcC=cGcC, PQSfinder = PQSfinder, RunComposition = RunComposition)
          if(nrow(PQSM3a)>0){
            PQSM3a$Strand <- "+"}
          if(nrow(PQSM3aCOMP)>0){
            PQSM3aCOMP$Strand <- "-"}
          if(nrow(PQSM3aCOMP)>0 & nrow(PQSM3a)>0){
            PQSM3a <- rbind(PQSM3a, PQSM3aCOMP)}
          if(nrow(PQSM3aCOMP)>0 & nrow(PQSM3a) == 0){
            PQSM3a <- PQSM3aCOMP}
          if(nrow(PQSM3aCOMP)== 0 & nrow(PQSM3a) > 0){
            PQSM3a <- PQSM3a}
          if(nrow(PQSM3aCOMP)== 0 & nrow(PQSM3a) == 0){
            PQSM3a <- data.frame()}
          rm(PQSM3aCOMP)
        }
        if(nrow(PQSM3a)>0){
          PQSM3a <- PQSM3a[order(PQSM3a$Start, decreasing = FALSE),]
          row.names(PQSM3a) <- seq(1:nrow(PQSM3a))
        }})[3])
      if(Verborrea == TRUE){cat(
        paste0(" & Ended. ", as.numeric(nrow(PQSM3a)), " results found"))
        }
      {
        if(nrow(PQSM3a)>0){
          if(Verborrea == TRUE){cat("\nM3 - Scoring: ")}
          ##Scoring system
          #G4hunter calculus
          if (G4hunter == TRUE){
            if(Verborrea == TRUE){cat("G4Hunter (I")}
            Process[length(Process)+1] <- "M3-G4Hunter"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3a <- .G4Hunter(df = PQSM3a, NCores = NCores)
            )[3])
            if(Verborrea == TRUE){cat("/O). ")}}
          #PQSfinder calculus
          if (PQSfinder == TRUE){
            if(Verborrea == TRUE){cat("PQSfinder (I")}
            Process[length(Process)+1] <- "M3-PQSfinder"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3a <- .PQSfinderfun(PQSM3a, RunComposition, Bt=Bt, Pb=Pb, Fm=Fm, Em=Em, Ts=Ts, Et=Et, Is=Is, Ei=Ei, Ls=Ls, ET = ET, MinNRuns= MinNRuns)
            )[3])
            if(Verborrea == TRUE){cat("/O). ")}
          }
          #cGcC
          if(cGcC == TRUE){
            if(Verborrea == TRUE){cat("cGcC (I")}
            Process[length(Process)+1] <- "M3-cGcC"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3a <- .cGcCfun(df = PQSM3a, RunComposition = RunComposition)
            )[3])
            if(Verborrea == TRUE){cat("/O). ")}}
          #mean Score
          if((G4hunter == TRUE & PQSfinder == TRUE)|(G4hunter == TRUE & cGcC == TRUE)|(PQSfinder == TRUE & cGcC == TRUE)){
            if(Verborrea == TRUE){ cat("MeanScore (I")}
            Process[length(Process)+1] <- "M3-meanScore"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3a <- .ScoreFun(df= PQSM3a, G4hunter = G4hunter, PQSfinder = PQSfinder, cGcC = cGcC, WeightParameters = WeightParameters)
            )[3])
            if(Verborrea == TRUE){cat("/O). ")}}

          if(Verborrea == TRUE){cat("\nQuantification Started")}
          ##Quantification
          {    Process[length(Process)+1] <- "M3-Quantification"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3a <- .LoopQuantification(df=PQSM3a, LoopSeq = LoopSeq))[3])
            if(Verborrea == TRUE){cat(" & Ended.")}}
          ##Qualification
          # Known G4 Sequences
          {
            if (((KnownQuadruplex== TRUE|KnownNOTQuadruplex == TRUE)|(KnownQuadruplex== TRUE & KnownNOTQuadruplex == TRUE)) & RunComposition == "G"|RunComposition == "C"){
              if(Verborrea == TRUE){cat("\nM3 - KnownQuadruplexStarted")}
              Process[length(Process)+1] <- "M3-KnownQuadruplex"
              Time[length(Time)+1] <- as.numeric(system.time(
                PQSM3a <- .KnownQuadFun(df = PQSM3a, DNA = DNA, KnownNOTQuadruplex = KnownNOTQuadruplex, KnownQuadruplex = KnownQuadruplex, RunComposition = RunComposition)
              )[3])
              if(Verborrea == TRUE){cat(" & Ended.")}}
            if(Verborrea == TRUE){cat("\nMethod 3B Started")}
            ## Method 3B by Frequency
            Process[length(Process)+1] <- "M3-M3B"
            Time[length(Time)+1] <- as.numeric(system.time(
              PQSM3b <- .M2B(df = PQSM3a, FreqWeight = FreqWeight, RunComposition = RunComposition)
            )[3])
            if(Verborrea == TRUE){cat(" & Ended.")}}
        }else{cat("\nNo results found with M3.")}
      }
    }
    ######  Times
    {
      FunTime <- data.frame(Process, Time, stringsAsFactors = FALSE)
      colnames(FunTime) <- c("Process", "sec.Time")
      FunTime$Date <- NA
      FunTime$Date[1] <- as.character(Date)
      FunTime$sec.Time[2] <-  sum(FunTime$sec.Time[3:nrow(FunTime)])
      FunTime$Date[2] <- as.character(Sys.time())
      if(FunTime$sec.Time[2] > 60) {
        FunTime$min.Time <-round(FunTime$sec.Time/60, 2)}
      if(FunTime$sec.Time[2] > 3600) {
        FunTime$hour.Time <-round(FunTime$sec.Time/3600,2)}
      if(FunTime$sec.Time[2] > 86400) {
        FunTime$day.Time <-round(FunTime$sec.Time/86400,2)}
    }
    ######  Config
    {        Variable <- c("Name","SeqSize" , "DNA", "Complementary", "RunComposition", "MaxRunSize", "MinRunSize",
                           "MaxNRuns", "MinNRuns", "BulgeSize", "MaxIL", "MaxPQSSize", "MinPQSSize", "MaxLoopSize",
                           "MinLoopSize", "LoopSeq", "Method2", "Method3", "G4hunter", "PQSfinder", "PQSfinder.Bt",
                           "PQSfinder.Pb", "PQSfinder.Fm", "PQSfinder.Em","PQSfinder.Ts", "PQSfinder.Is", "PQSfinder.Ls", "PQSfinder.Et",
                           "PQSfinder.Ei", "PQSfinder.ET", "cGcC", "WeightParameters", "FreqWeight", "KnownQuadruplex", "KnownNOTQuadruplex", "NCores", "SeqG%", "SeqC%",
                           "G4iMGrinder.Version", "GiG.DB.BioInformatic.Version", "GiG.DB.Refs.Version", "GiG.DB.BioPhysical.Version",
                           "M1A.Results", "M1B.Results")
      Value <- c(
        Name, ifelse(Complementary == TRUE, yes = nchar(Sequence)*2, no = nchar(Sequence)),
        DNA, Complementary, RunComposition, MaxRunSize, MinRunSize,
        MaxNRuns, MinNRuns, BulgeSize, MaxIL, MaxPQSSize, MinPQSSize, MaxLoopSize,
        MinLoopSize, paste0(LoopSeq, collapse = ", "), Method2, Method3, G4hunter, PQSfinder, Bt,
        Pb, Fm, Em, Ts, Is, Ls, Et, Ei,ET,  cGcC, paste0(WeightParameters, collapse = ", "), FreqWeight, KnownQuadruplex, KnownNOTQuadruplex, NCores,
        round(100*ifelse(Complementary == TRUE,
                         no =  str_count(Sequence,pattern = "G")/nchar(Sequence),
                         yes =  (str_count(Sequence,pattern = "C") + str_count(Sequence,pattern = "G"))/(nchar(Sequence)*2)), 1),
        round(100*ifelse(Complementary == TRUE,
                         no =  str_count(Sequence,pattern = "C")/nchar(Sequence),
                         yes =  (str_count(Sequence,pattern = "C") + str_count(Sequence,pattern = "G"))/(nchar(Sequence)*2)), 1),
        packageDescription("G4iMGrinder")$Version,
        GiG.DB$Version$Version[1:3],
        as.numeric(nrow(PQSM1a) + ifelse(Complementary == TRUE, yes = nrow(PQSM1aCOMP), no = 0)),
        as.numeric(nrow(PQSM1b) + ifelse(Complementary == TRUE, yes = nrow(PQSM1bCOMP), no = 0))
        )
      Configuration <- data.frame(Variable)
      Configuration$Value <- Value

    }
    ######  Finishing up
    {
      PQSM1a <- NULL
      PQSM1b <- NULL
      PQSM1aCOMP <- NULL
      PQSM1bCOMP <- NULL

      dfs <- Filter(function(x) is(x, "data.frame"), mget(ls()))
    }
    if(Verborrea == TRUE){
    if(RunComposition == "C"){
      cat("\nRemember, i-Motifs are puntuated inversily to normal PQS. -100 scores are the best most probable to form i-Motifs.")}
    cat(paste0("\nAnalysis of ", Name, " finished in ", round(FunTime$sec.Time[2],0)," s.\n Have a great day :) EBR\n"))
    beepr::beep(4)}
    return(dfs)
  } else {
    cat("Please Select a method to analyze the data (Method 2 and/or Method 3).")}
}
################################################################################
#########  G4-iM Grinder: End of function.
################################################################################
####
####
################################################################################
#########  GiGList.Analysis: Analysis of G4-iM Grinder Results
################################################################################
########
GiGList.Analysis <- function(GiGList, iden, Density = 100000, ScoreMin = 40, FreqMin = 10,  LengthMin = 50,  Verborrea = TRUE){

  # Data Table
  ResultTable <- data.frame(matrix(vector(), 0, 1,
                                   dimnames=list(c(), c("Name"))),
                            stringsAsFactors=F)

  # Search Info
  require(stringr)
  Name <- as.character(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Name"])
  SeqSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="SeqSize"])
  DNA <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="DNA"])
  Complementary <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Complementary"])
  Method2 <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Method2"])
  Method3 <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Method3"])
  MinPQSSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MinPQSSize"])
  G <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="SeqG%"])
  C <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="SeqC%"])
  RunComposition <- as.character(GiGList$Configuration$Value[GiGList$Configuration$Variable =="RunComposition"])

  # Sequence part
  {
    ResultTable[length(ResultTable$Name) +1,] <- NA
    ResultTable$Name[length(ResultTable$Name)] <- Name
    ResultTable$iden[length(ResultTable$Name)] <- iden
    LengthTotal <- SeqSize
    ResultTable$Length[length(ResultTable$Name)] <- LengthTotal
    ResultTable$SeqG <- G
    ResultTable$SeqC <- C
  }

  #  Results with densities
  if(Method2 == TRUE){
    ResultTable$nM2a[length(ResultTable$Name)] <- as.numeric(nrow(GiGList$PQSM2a))
    ResultTable$nM2a.D[length(ResultTable$Name)] <- round((ResultTable$nM2a[length(ResultTable$Name)]/LengthTotal)*Density,2)
    ResultTable$nM2b[length(ResultTable$Name)] <- as.numeric(nrow(GiGList$PQSM2b))
    ResultTable$nM2b.D[length(ResultTable$Name)] <- round((ResultTable$nM2b[length(ResultTable$Name)]/LengthTotal)*Density,2)
    if(nrow(GiGList$PQSM2a) == 0){
      ResultTable$nM2b[length(ResultTable$Name)] <- 0
      ResultTable$nM2b.D[length(ResultTable$Name)] <-  0
      ResultTable$nM2a.S[length(ResultTable$Name)] <- 0
      ResultTable$nM2a.D.S[length(ResultTable$Name)]  <- 0
      ResultTable$nM2b.S[length(ResultTable$Name)] <- 0
      ResultTable$nM2b.D.S[length(ResultTable$Name)]  <- 0
      ResultTable$nM2b.F[length(ResultTable$Name)] <- 0
    } else {
      if(RunComposition == "C") {
        Index <- GiGList$PQSM2a$Score <= ScoreMin
        ResultTable$nM2a.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM2a.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM2a.S[length(ResultTable$Name)]/LengthTotal)*Density,2)
        Index <- GiGList$PQSM2b$Score <= ScoreMin
        ResultTable$nM2b.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM2b.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM2b.S[length(ResultTable$Name)]/LengthTotal)*Density,2)

        Index <- GiGList$PQSM2b$freq >= FreqMin
        ResultTable$nM2b.F[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
      } else {
        Index <- GiGList$PQSM2a$Score >= ScoreMin
        ResultTable$nM2a.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM2a.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM2a.S[length(ResultTable$Name)]/LengthTotal)*Density,2)
        Index <- GiGList$PQSM2b$Score >= ScoreMin
        ResultTable$nM2b.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM2b.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM2b.S[length(ResultTable$Name)]/LengthTotal)*Density,2)

        Index <- GiGList$PQSM2b$freq >= FreqMin
        ResultTable$nM2b.F[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
      }
    }
  }
  if(Method3 == TRUE){
    ResultTable$nM3a[length(ResultTable$Name)] <- as.numeric(nrow(GiGList$PQSM3a))
    ResultTable$nM3a.D[length(ResultTable$Name)] <- round((ResultTable$nM3a[length(ResultTable$Name)]/LengthTotal)*Density,2)
    ResultTable$nM3b[length(ResultTable$Name)] <- as.numeric(nrow(GiGList$PQSM3b))
    ResultTable$nM3b.D[length(ResultTable$Name)] <- round((ResultTable$nM3b[length(ResultTable$Name)]/LengthTotal)*Density,2)

    if(nrow(GiGList$PQSM3a) == 0){
      ResultTable$nM3b[length(ResultTable$Name)] <- 0
      ResultTable$nM3b.D[length(ResultTable$Name)] <- 0
      ResultTable$nM3a.S[length(ResultTable$Name)] <- 0
      ResultTable$nM3a.D.S[length(ResultTable$Name)]  <- 0
      ResultTable$nM3b.S[length(ResultTable$Name)] <- 0
      ResultTable$nM3b.D.S[length(ResultTable$Name)]  <- 0
      ResultTable$nM3b.F[length(ResultTable$Name)] <- 0
      ResultTable$nM3b.L[length(ResultTable$Name)] <- 0
    } else {

      if(RunComposition == "C") {
        Index <- GiGList$PQSM3a$Score <= ScoreMin
        ResultTable$nM3a.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM3a.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM3a.S[length(ResultTable$Name)]/LengthTotal)*Density,2)
        Index <- GiGList$PQSM3b$Score <= ScoreMin
        ResultTable$nM3b.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM3b.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM3b.S[length(ResultTable$Name)]/LengthTotal)*Density,2)

        Index <- GiGList$PQSM3b$freq >= FreqMin
        ResultTable$nM3b.F[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        Index <- GiGList$PQSM3a$Length >= LengthMin
        ResultTable$nM3a.L[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
      } else {
        Index <- GiGList$PQSM3a$Score >= ScoreMin
        ResultTable$nM3a.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM3a.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM3a.S[length(ResultTable$Name)]/LengthTotal)*Density,2)
        Index <- GiGList$PQSM3b$Score >= ScoreMin
        ResultTable$nM3b.S[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        ResultTable$nM3b.D.S[length(ResultTable$Name)]  <- round((ResultTable$nM3b.S[length(ResultTable$Name)]/LengthTotal)*Density,2)

        Index <- GiGList$PQSM3b$freq >= FreqMin
        ResultTable$nM3b.F[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
        Index <- GiGList$PQSM3a$Length >= LengthMin
        ResultTable$nM3a.L[length(ResultTable$Name)] <- as.numeric(length(Index[Index == TRUE]))
      }

    }
  }

  ##Config protocols
  {
    ResultTable$Config[length(ResultTable$Name)] <- paste0(
      paste0("|", "(D)Density(", Density, ")|"),
      if(!is.na(ScoreMin)){ if (is.numeric(ScoreMin)){ if(ScoreMin != ""){
        paste0("(S)Score=>", ScoreMin, "|")}else{"."} }else{"."} }else{"."},
      if (!is.na(FreqMin)){ if (is.numeric(FreqMin)){if(FreqMin != ""){
        paste0("(F)Freq.=>", FreqMin, "|")}else{"."} }else{"."} }else{"."},
      if (!is.na(LengthMin)){ if (is.numeric(LengthMin)){if(LengthMin != ""){if(LengthMin > MinPQSSize){
        paste0("(L)Length=>", LengthMin, "|\n")}else{"."}}}})
  }
  if(Verborrea == TRUE){
    cat("\nG4-iM Grinder - Description of the results summary table:")
    cat("\nThe first 7 columns show the sequence size and composition")
    cat("\nColumns nM2 and nM3 are the n? of structures found with Methods 2 & 3")
    cat(paste0("\nD is the total structure Density per ", Density, " nucleotides of the sequence"))
    if(!is.na(ScoreMin)){if (is.numeric(ScoreMin)){if(ScoreMin != ""){
      cat(paste0("\nS is the total structures found with Score > ", ScoreMin))}}}
    if (!is.na(FreqMin)){ if (is.numeric(FreqMin)){if(FreqMin != ""){
      cat(paste0("\nF is the total structures found with frequency >", FreqMin))}}}
    if (!is.na(LengthMin)){ if (is.numeric(LengthMin)){if(LengthMin != ""){if(LengthMin > MinPQSSize){
      cat(paste0("\nL is the total structures found with the Length >", LengthMin, " nucleotides\n"))
    }}}}
  }
  return(ResultTable)
}
########
################################################################################
#########  GiG.M3Structure: Analysis of M3A candidates - Higher Order Structure potential subunits
################################################################################
########
GiG.M3Structure <- function(GiGList, M3ACandidate, MAXite){
  #Preparing DATA
  {
    #M3A result


    BBB <- GiGList$PQSM3a[M3ACandidate,]
    stopifnot(nrow(BBB) > 0 & nrow(BBB) < 2)

    #M2A result
    if("Chromosome" %in% colnames(BBB)){
      AAA <- GiGList$PQSM2a[
        GiGList$PQSM2a$Chromosome == BBB$Chromosome &
          GiGList$PQSM2a$Start >= BBB$Start &
          GiGList$PQSM2a$Start <=  BBB$Finish
        ,]
    } else {
      AAA <- GiGList$PQSM2a[
        GiGList$PQSM2a$Start >= BBB$Start & GiGList$PQSM2a$Start <=  BBB$Finish ,]
    }
    if("Strand" %in% colnames(BBB)){
      AAA <- AAA[AAA$Strand ==BBB$Strand,]
    }

    stopifnot(nrow(AAA) > 0)
    AAA$RES <- row.names(AAA)
    row.names(AAA) <- seq(1, nrow(AAA), 1)
  }

  #selecting invalids that overlap, FUNCTION
  Overlapfun <- function(AAA, i){
    index2 <- (
      (   AAA$Start [i] >= AAA$Start  &
            AAA$Start [i] <= AAA$Finish &
            AAA$Finish[i] >= AAA$Start  &
            AAA$Finish[i] >= AAA$Finish)  |

        (   AAA$Start [i] <= AAA$Start  &
              AAA$Start [i] <= AAA$Finish &
              AAA$Finish[i] >= AAA$Start  &
              AAA$Finish[i] >= AAA$Finish)  |

        (     AAA$Start [i] <= AAA$Start  &
                AAA$Start [i] <= AAA$Finish &
                AAA$Finish[i] >= AAA$Start  &
                AAA$Finish[i] <= AAA$Finish) |

        (AAA$Start [i] >= AAA$Start  &
             AAA$Start [i] <= AAA$Finish &
             AAA$Finish[i] >= AAA$Start  &
            AAA$Finish[i] <= AAA$Finish)

    )
    AAA$OK[index2] <- "X"
    AAA$OK[i]  <- "PQS"
    return(AAA)}

  #Highest score sequencial analysis function
  HSAfun <- function(AAA){
    HSA <- AAA
    HSA$SxS <- round(x = HSA$Score*ifelse(HSA$Conf.Quad.Seqs == "", 1, 1.5),1)
    HSA$OK <- NA

    repeat {
      HSA2 <- HSA[is.na(HSA$OK),]
      i <- as.numeric(row.names(HSA2[HSA2$SxS == max(HSA2$SxS),]))
      if(length(i) >1){
            i <- i[sample(x = 1:length(i), size = 1)]
      }
      HSA <- Overlapfun(AAA = HSA, i = i)
      if(nrow(HSA[is.na(HSA$OK),]) < 1){
        break()
      }
    }
    HSA <- HSA[HSA$OK =="PQS",]
    return(HSA) }
  HSA1 <- HSAfun(AAA)

  # Random structure conformation function
  RSCfun <- function(AAA, MAXite){

    AAA$SxS <- 0
    AAA$OK <- NA
    Options <- data.frame(Ite = numeric(0),nPQS = numeric(0), MeanScore = numeric(0), PQSLengthPercent= numeric(0))
    ite <- 1
    Options[ite,] <- NA
    Options$idenPQS <- NA
    repeat{
      RSC <- AAA
      Options[ite,] <- NA

      repeat {
        RSC2 <- RSC[is.na(RSC$OK),]
        i <- as.numeric(row.names(RSC2[sample(x = nrow(RSC2),size = 1),]))
        RSC <- Overlapfun(AAA = RSC, i = i)
        rm(RSC2)
        if(nrow(RSC[is.na(RSC$OK),]) < 1){
          break()
        }

      }

      RSC <- RSC[RSC$OK == "PQS",]

      Options$Ite[ite] <- ite
      Options$MeanScore[ite] <- round(mean(RSC$Score),1)
      Options$nPQS[ite] <- nrow(RSC)
      Options$PQSLengthPercent[ite] <- round(100*sum(as.numeric(RSC$Length))/BBB$Length,1)
      Options$idenPQS[ite] <- as.character(paste0(row.names(RSC), collapse = " "))
      ite <- ite +1
      if(ite > MAXite){
        break()}
    }
    Options$nMS <- round(Options$MeanScore*Options$PQSLengthPercent/100,1)

    return(Options)

  }
  RSC <- RSCfun(AAA, MAXite = MAXite)

  #unique Random conformation and extraction of RAH and RAnH
  RSC$Ite <- NULL
  uRSC <- unique(RSC)
  Arrangements <- NULL
  Arrangements$HSA <- uRSC[uRSC$idenPQS == paste0(row.names(HSA1), collapse = " "),]
  Arrangements$RAH <-uRSC[uRSC$MeanScore == max(uRSC$MeanScore),]
  Arrangements$RAnH <- uRSC[uRSC$nMS == max(uRSC$nMS),]


  Ending <- list(AAA, BBB, uRSC, Arrangements)
  names(Ending) <- c("M2", "M3", "Potential.Arrangements", "Best.Arrangements")
  return(Ending)
}
########
################################################################################
#########  GiGList.Updater: Updater function for an existant G4-iM Grinder Result
################################################################################
########
GiGList.Updater <- function(GiGList,  ChangeRunComposition = F, LoopSeq= c("G", "T", "A", "C"),
                           WeightParameters=c(50,50,0), G4hunter=T,  PQSfinder=T, Bt=14, Pb=17, Fm=3, Em=1, Ts=4, Et=1, Is=-19, Ei=1, Ls=-16, ET = 1, FreqWeight=0,
                           KnownQuadruplex= T, KnownNOTQuadruplex= T,  RunFormula = F, NCores = 1){
  ##Loading packages
  .PackageLoading()

  GiGList$Configuration$Variable <- as.character(GiGList$Configuration$Variable)
  Name <- as.character(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Name"])
  SeqSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="SeqSize"])
  DNA <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="DNA"])
  Complementary <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Complementary"])
  RunComposition <- as.character(GiGList$Configuration$Value[GiGList$Configuration$Variable =="RunComposition"])
  MaxRunSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MaxRunSize"])
  MinRunSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MinRunSize"])
  MaxNRuns <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MaxNRuns"])
  MinNRuns <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MinNRuns"])
  BulgeSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="BulgeSize"])
  MaxIL <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MaxIL"])
  MaxPQSSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MaxPQSSize"])
  MinPQSSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MinPQSSize"])
  MaxLoopSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MaxLoopSize"])
  MinLoopSize <- as.numeric(GiGList$Configuration$Value[GiGList$Configuration$Variable =="MinLoopSize"])
  Method2 <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Method2"])
  Method3 <- as.logical(GiGList$Configuration$Value[GiGList$Configuration$Variable =="Method3"])
  cGcC <- F

  #Inverting RunComposition
  if(ChangeRunComposition == T & (RunComposition == "G"| RunComposition == "C")){
    if(Complementary == TRUE){

      RunComposition <- ifelse(test = RunComposition == "G", yes =  "C", no =  "G")

    SeqSize <- SeqSize/2
    if(Method2 == TRUE){
      GiGList$PQSM2a$Sequence <- G4iMGrinder:::.CompSeqFun(DNA = DNA, Sequence = GiGList$PQSM2a$Sequence)
      StartNew <- SeqSize - GiGList$PQSM2a$Finish+1
      FinishNew <- SeqSize - GiGList$PQSM2a$Start+1
      GiGList$PQSM2a$Start <- StartNew
      GiGList$PQSM2a$Finish <- FinishNew
      GiGList$PQSM2a$Strand <- ifelse(test = GiGList$PQSM2a$Strand == "+", yes =  "-", no =  "+")
      GiGList$PQSM2a <- select(GiGList$PQSM2a, Start, Finish, Length, Runs, IL, mRun, Sequence,  Strand )
      GiGList$PQSM2a <- arrange(GiGList$PQSM2a, Start)
      rownames(GiGList$PQSM2a) <- seq(1, nrow(GiGList$PQSM2a), 1)
    }
    if(Method3 == TRUE){
      GiGList$PQSM3a$Sequence <- G4iMGrinder:::.CompSeqFun(DNA = DNA, Sequence = GiGList$PQSM3a$Sequence)
      StartNew <- SeqSize - GiGList$PQSM3a$Finish+1
      FinishNew <- SeqSize - GiGList$PQSM3a$Start+1
      GiGList$PQSM3a$Start <- StartNew
      GiGList$PQSM3a$Finish <- FinishNew
      GiGList$PQSM3a$Strand <- ifelse(test = GiGList$PQSM3a$Strand == "+", yes =  "-", no =  "+")
      GiGList$PQSM3a <- select(GiGList$PQSM3a, Start, Finish, Length, Runs, IL, mRun, Sequence, Strand )
      GiGList$PQSM3a <- arrange(GiGList$PQSM3a, Start)
      rownames(GiGList$PQSM3a) <- seq(1, nrow(GiGList$PQSM3a), 1)
    }
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="RunComposition"] <- RunComposition
    } else {
      cat("\n RunComposition inversion cannot be done on a single genomic strand. The genomic complementary analysis is required.")
    }
  }

  ## Calls to apply
  if (Method2 == TRUE){
    # Score
    if(G4hunter+PQSfinder>0){
      if(G4hunter == TRUE){
        GiGList$PQSM2a <- G4iMGrinder:::.G4Hunter(df = GiGList$PQSM2a, NCores = NCores)
      }
      if(PQSfinder == TRUE){
        GiGList$PQSM2a <- G4iMGrinder:::.PQSfinderfun(GiGList$PQSM2a, RunComposition= RunComposition, Bt=Bt, Pb=Pb, Fm=Fm, Em=Em, Ts=Ts, Et=Et, Is=Is, Ei=Ei, Ls=Ls, ET = ET, MinNRuns= MinNRuns)
      }
      GiGList$PQSM2a <- G4iMGrinder:::.ScoreFun(df= GiGList$PQSM2a, G4hunter = G4hunter, PQSfinder = PQSfinder, cGcC = F, WeightParameters = WeightParameters)
    }
    #Loop Quantification
    if(length(LoopSeq) >0){
      GiGList$PQSM2a <- G4iMGrinder:::.LoopQuantification(df=GiGList$PQSM2a, LoopSeq = LoopSeq)
    }
    # KnownQuadFun
    if(KnownQuadruplex + KnownNOTQuadruplex > 0){
      GiGList$PQSM2a <- .KnownQuadFun(GiGList$PQSM2a, KnownNOTQuadruplex = KnownNOTQuadruplex , KnownQuadruplex = KnownQuadruplex, RunComposition = RunComposition, DNA = DNA)
    }
    GiGList$PQSM2b <- G4iMGrinder:::.M2B(df = GiGList$PQSM2a, FreqWeight = FreqWeight, RunComposition = RunComposition)
  }
  if (Method3 == TRUE){
    # Score
    if(G4hunter +PQSfinder > 0){
      if(G4hunter == TRUE){
        GiGList$PQSM3a <- .G4Hunter(df = GiGList$PQSM3a, NCores = NCores)
      }
      if(PQSfinder == TRUE){
        GiGList$PQSM3a <- .PQSfinderfun(GiGList$PQSM3a, RunComposition, Bt=Bt, Pb=Pb, Fm=Fm, Em=Em, Ts=Ts, Et=Et, Is=Is, Ei=Ei, Ls=Ls, ET = ET, MinNRuns= MinNRuns)
      }
      GiGList$PQSM3a <- .ScoreFun(df= GiGList$PQSM3a, G4hunter = G4hunter, PQSfinder = PQSfinder, cGcC = F, WeightParameters = WeightParameters)
    }
    #Loop Quantification
    if(length(LoopSeq) >0){
      GiGList$PQSM3a <- .LoopQuantification(df=GiGList$PQSM3a, LoopSeq = LoopSeq)
    }
    # KnownQuadFun
    if(KnownQuadruplex + KnownNOTQuadruplex >0){
      GiGList$PQSM3a <- .KnownQuadFun(GiGList$PQSM3a, KnownNOTQuadruplex = KnownNOTQuadruplex , KnownQuadruplex = KnownQuadruplex, RunComposition = RunComposition, DNA = DNA)
    }
    GiGList$PQSM3b <- .M2B(df = GiGList$PQSM3a, FreqWeight = FreqWeight, RunComposition = RunComposition)
  }

  ##Changing variables of configuration
  if("Score" %in% colnames(GiGList$PQSM2a)){
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="WeightParameters"] <- paste0(WeightParameters, collapse = " ")
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="FreqWeight"] <- FreqWeight
  }
  if("G4Hunter" %in% colnames(GiGList$PQSM2a)){
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="G4hunter"] <- TRUE
  }
  if("pqsfinder" %in% colnames(GiGList$PQSM2a)){
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Bt"] <- Bt
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Pb"] <- Pb
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Fm"] <- Fm
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Em"] <- Em
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Ts"] <- Ts
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Et"] <- Et
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Is"] <- Is
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Ei"] <- Ei
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.Ls"] <- Ls
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="PQSfinder.ET"] <- ET
    }
  if("Conf.Quad.Seqs" %in% colnames(GiGList$PQSM2a)){
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="KnownQuadruplex"] <- TRUE
    }
  if("Conf.NOT.Quad.Seqs" %in% colnames(GiGList$PQSM2a)){
    GiGList$Configuration$Value[GiGList$Configuration$Variable =="KnownNOTQuadruplex"] <- TRUE
  }

  return(GiGList)
  }
########
########
########
########
########
########
########
########
########
########
########
########
########
########
########
########
