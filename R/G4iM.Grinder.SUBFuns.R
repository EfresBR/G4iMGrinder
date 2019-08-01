
### SUBFUNCTIONS OF G4iM GRINDER

{
  # Sequence
  # Complementary calculations
  .CompSeqFun <- function(DNA, Sequence){
    ifelse(DNA == TRUE, yes = Sequence2 <- stringi::stri_reverse(chartr(old = "ATCG", new = "TAGC", x = Sequence)),
           no = Sequence2 <- stringi::stri_reverse(chartr(old = "AUCG", new = "UAGC", x = Sequence)))
    return(Sequence2)
  }
  #
  ##### Method 1
  # Search for Runs. At first it will search for perfect runs within the Min and Max RunSize. No overlapping within these runs is accepted.
  # Hence GGGG will only be that and not GGG as well.
  # If BulgeSize >0, then it will seek to find runs which fit Min and Max RunSize with BulgeSize. However, for this it will not take
  # into account the perfect runs previously found. This is to avoid more stable versions of Gruns with less Loops, besides improving performance.
  # Hence, a sequence which contains GTGGG, first will find GGG and then it will ?NOT be associated to the previous G with the T loop.
  # However, GTGG, will be identified as a Run with 1 Loop.
  .M1A<- function(RunComposition, MaxRunSize, MinRunSize , BulgeSize, Sequence){
    stopifnot(nchar(Sequence)>0)
    pattern <- NULL
    Glocations <- NULL
    Start <- NULL
    Finish <- NULL
    SequenceALT <- Sequence

    ### Perfect Runs
    {
      for (i in MaxRunSize:MinRunSize){
        pattern <- str_c(rep(RunComposition, i), collapse = "")
        Glocations <- gregexpr(pattern, text = SequenceALT)
        Start <- c(Start, as.numeric(Glocations[[1]]))
        Finish <- c(Finish, as.numeric(Glocations[[1]]) + i-1)
        SequenceALT <- str_replace_all(SequenceALT, pattern = pattern, replacement = str_c(rep("*", i), collapse = ""))
        # modyfying as to make faster for IL search, and smaller tetrad searchs non overlapping.
      }
      Score <- data.frame(Start, Finish)
      Index <- Score$Start >0
      Score <- Score[Index,]
      if(nrow(Score) >0){
        Score$IL <- 0}
    }
    ### Imperfect Runs
    if(BulgeSize > 0){
      Glocations <- str_locate_all(pattern = RunComposition, string = SequenceALT )
      Glocations <- Glocations[[1]]
      Glocations <- Glocations[,1]
      if(length(Glocations) >= MinRunSize){
        count <- 0
        repeat {
          n <- MinRunSize + count
          df <- data.frame(Start = Glocations, Finish = NA, IL = NA)
          index <- Glocations[c((n):length(Glocations), rep(NA, n-1))]-Glocations < (n+BulgeSize)
          index[is.na(index)] <- F
          if(sum(index) == 0){
            break
          }
          df$Finish <- df$Start[c((n):nrow(df), rep(NA, n-1))]
          df <- df[index,]
          df$IL <- df$Finish-df$Start-n+1
          Score <- rbind(Score, df)
          count <- count +1
          if(MinRunSize + count > MaxRunSize){
            break
          }
          if(n == length(Glocations)){
            break
          }
        }
      }
    }

    ### Safeguarding
    if(nrow(Score) > 0){
      Score <- unique(Score)
      Score <- Score[order(Score$Start, decreasing = FALSE),]
      rownames(Score) <- seq(1:nrow(Score))
      return(Score)
    } else {
      Score <- data.frame(Start = numeric(), Finish =numeric(), IL = numeric())
      return(Score)
    }

  }
  #
  # Will find the assosiation between the Runs identified previously. and Link each run to another Run when Min and Max LoopSizes are respected.
  # It will delete Runs which are not associated with other Runs.
  .M1B <- function(df, MinLoopSize, MaxLoopSize, MinRunSize, MinNRuns){
    if(nrow(df)>= MinNRuns){
      PQSM1b <- df

      # Loops between Gruns.
      PQSM1b$Pre <- NA
      PQSM1b$Pos <- NA
      PQSM1b$AAA <- c(PQSM1b$Start[2:nrow(PQSM1b)], PQSM1b$Start[nrow(PQSM1b)])
      PQSM1b$Pos <-  PQSM1b$AAA -PQSM1b$Finish -1
      PQSM1b$Pos[nrow(PQSM1b)] <- MaxLoopSize*10
      PQSM1b$AAA <- NULL
      PQSM1b$Pre <- c(10*MaxLoopSize, PQSM1b$Pos[1:nrow(PQSM1b)-1])
      Index <- (PQSM1b$Pre > MaxLoopSize) & (PQSM1b$Pos > MaxLoopSize)
      PQSM1b <- PQSM1b[!Index,]

      # Trimming GrunSize to fit MinLoopSize and MinRunSize, by POSterior loop only for non IL- runs
      PQSM1b$RunSize <- PQSM1b$Finish-PQSM1b$Start +1
      A <- NULL

      # Checking length of run to try and reduce it length (respecting MinRunSize) to cede space to minLoop in Pos.
      Index <-(PQSM1b$Pos < MinLoopSize) & (PQSM1b$IL == 0 )

      # Fixing those that fit
      A <- PQSM1b$RunSize[Index] - (MinLoopSize-PQSM1b$Pos[Index]) - (MinRunSize)
      AIndex <- (A > 0)
      PQSM1b$Finish[Index][AIndex]  <- PQSM1b$Finish[Index][AIndex]-A[AIndex]
      PQSM1b$Pos[Index][AIndex]  <- PQSM1b$Pos[Index][AIndex]+A[AIndex]

      # These do not fit MinLoop but can get reduced, so they will be in the hope of Pre modifictation to fit.
      BIndex <- (PQSM1b$RunSize[Index][!AIndex]  >  MinRunSize)
      PQSM1b$Pos[Index][!AIndex][BIndex]<- PQSM1b$Pos[Index][!AIndex][BIndex] + PQSM1b$RunSize[Index][!AIndex][BIndex]-MinRunSize
      PQSM1b$RunSize[Index][!AIndex][BIndex] <- MinRunSize
      PQSM1b$Finish[Index][!AIndex][BIndex] <- PQSM1b$Start[Index][!AIndex][BIndex] + (MinRunSize-1)
      PQSM1b$Pre <- c(2*MaxLoopSize, PQSM1b$Pos[1:nrow(PQSM1b)-1])

      # Checking Pre to try and reduce Run size to cede space to minLoop.
      Index <-(PQSM1b$Pre < MinLoopSize) & (PQSM1b$IL == 0 )
      A <- PQSM1b$RunSize[Index] - (MinLoopSize-PQSM1b$Pre[Index]) - (MinRunSize)
      AIndex <- (A > 0)
      PQSM1b$Start[Index][AIndex]  <- PQSM1b$Start[Index][AIndex]+A[AIndex]
      PQSM1b$Pre[Index][AIndex]  <- PQSM1b$Pos[Index][AIndex]+A[AIndex]
      PQSM1b$RunSize <- PQSM1b$Finish-PQSM1b$Start +1
      PQSM1b$Pos <- c(PQSM1b$Pre[2:nrow(PQSM1b)],2*MaxLoopSize)
      PQSM1b$RunSize <- NULL

      # Linking to column number the ones that can be linked to neighbors
      Index <- (PQSM1b$Pos >= MinLoopSize) & (PQSM1b$Pos <= MaxLoopSize)
      row.names(PQSM1b)<- seq(1:nrow(PQSM1b))
      PQSM1b$Link <- as.numeric(row.names(PQSM1b))+1
      PQSM1b$Link[!Index] <- NA
      PQSM1b$Link[nrow(PQSM1b)] <- NA
      PQSM1b$FK <- as.numeric(row.names(PQSM1b))

      #Checking if the unfit can get linked to other runs which are not direct neighbors
      if (nrow(PQSM1b[PQSM1b$Pos < MinLoopSize,]) > 0){
        PQSM1b$AAA <- c(PQSM1b$Start[2:nrow(PQSM1b)], PQSM1b$Start[nrow(PQSM1b)]+500)
        count <- 1
        while(nrow(PQSM1b[PQSM1b$Pos < MinLoopSize,]) > 0){
          count <- count +1
          PQSM1b$AAA <- c(PQSM1b$AAA[2:nrow(PQSM1b)], (PQSM1b$Finish[nrow(PQSM1b)]+500))
          PQSM1b$Pos[!Index] <- PQSM1b$AAA[!Index]-PQSM1b$Finish[!Index]-1
          Index2 <- (PQSM1b$Pos[!Index] >= MinLoopSize) & (PQSM1b$Pos[!Index] <= MaxLoopSize)
          PQSM1b$Link[!Index][Index2] <- (PQSM1b$FK[!Index][Index2] + count)

          Index <- (PQSM1b$Pos >= MinLoopSize) & (PQSM1b$Pos <= MaxLoopSize)
          if(count == 10) {break()}
        }
        PQSM1b$AAA <- NULL
      }
      return(PQSM1b)
    } else {
      PQSM1b <- data.frame(Start = numeric(), Finish =numeric(), IL = numeric(),
                           Pre = numeric(), Pos = numeric(), Link = numeric(), FK = numeric())
      return(PQSM1b)}
  }
  #
  ##### Method 2
  # Detects PQS (or i-motifs etc) which obey Max and MinPQSSize and Max and Min n? of Runs. MaxGrun size is not neccesary.
  .M2a <- function(df, RunComposition, MinPQSSize, MinNRuns, MaxPQSSize, MaxNRuns, PQSfinder, RunFormula = FALSE, Sequence, MaxIL, NCores, Verborrea, BulgeSize, G4hunter,  cGcC){
    if(nrow(df)>=MinNRuns){
      PQSM2 <- df

      #Search for distant G runss from maxPQSSize
      if(NCores == 1){
        Start <- NULL
        Finish <- NULL
        Runs <- NULL
        IL <- NULL
        mRun <- NULL
        mRuns <- NULL
        count <- NULL
        Length <- NULL
        IL <- NULL
        ILs <- NULL
        Formula <- NULL
        Formulas <- NULL
        Index <- PQSM2$Link >0
        for (i in 1:nrow(PQSM2[Index,])){
          count <- PQSM2$Link[i]
          if(is.na(count)){
            next()
          }
          Length <- 0
          Run <- 1
          if(BulgeSize >0){
            ILs <- PQSM2$IL[i]}
          if(PQSfinder == TRUE){
            mRuns <- PQSM2$Finish[i]- PQSM2$Start[i]+1 - ifelse(BulgeSize == 0, 0, PQSM2$IL[i])}
          if(RunFormula == TRUE){
            Formulas <- paste0(PQSM2$Start[i]+1-PQSM2$Start[i],":", PQSM2$Finish[i]+1-PQSM2$Start[i])}
          repeat {
            Length <- PQSM2$Finish[count]- PQSM2$Start[i]+1
            if(Length > MaxPQSSize){
              break()
            }
            Run <- Run +1
            if(MinNRuns <= MaxNRuns){
              if(Run > MaxNRuns) {break()}
            }
            if(BulgeSize >0){
              ILs<- ILs + PQSM2$IL[count]
              if(ILs > MaxIL)
                break()}
            if(PQSfinder == TRUE){
              mRuns <- mRuns + (PQSM2$Finish[count]- PQSM2$Start[count]+1 - ifelse(BulgeSize >0, PQSM2$IL[count], 0))}
            if(RunFormula == TRUE){
              Formulas <- paste0(Formulas,".", PQSM2$Start[count]+1-PQSM2$Start[i],":", PQSM2$Finish[count]+1-PQSM2$Start[i])}
            if(Run >= MinNRuns & Length >= MinPQSSize){
              Start[length(Start)+1] <- PQSM2$Start[i]
              Finish[length(Finish)+1] <-PQSM2$Finish[count]
              Runs[length(Runs)+1] <- Run
              if(BulgeSize >0){
                IL[length(IL)+1] <- ILs}
              if(PQSfinder == TRUE) {
                mRun[length(mRun)+1] <- mRuns/Run}
              if(RunFormula == TRUE){
                Formula[length(Formula)+1] <- Formulas}
            }
            count <- PQSM2$Link[count]
            if(is.na(count)== TRUE){
              break()
            }
          }
        }
        if(RunFormula == TRUE){
          PQSM2a <- as.data.frame(cbind(Start, Finish, Runs, IL, mRun,Formula))} else {
            PQSM2a <- as.data.frame(cbind(Start, Finish, Runs, IL, mRun))
          }
      } else {
        Index <- is.na(PQSM2$Link)
        Index <- PQSM2[!Index,]
        Index <- Index[complete.cases(Index),]

        M2Loop <- function(Line){
          Line <- as.numeric(Line)
          Start <- Line[1]
          Finish <- Line[2]
          IL <- Line[3]
          Link <- Line[6]

          V1 <- data.frame(Start = integer(), Finish = integer(), Runs = integer(), IL = integer(), mRun = integer(), Formula = character())

          count <- Link
          if(!is.na(count)){
            Length <- 0
            Run <- 1
            mRuns <- Finish- Start+1 - IL
            if(RunFormula == TRUE){
              Formulas <- paste0(1,":", Finish+1-Start)}

            repeat {
              Length <- PQSM2$Finish[count]- Start + 1
              if(Length > MaxPQSSize){
                break()
              }
              Run <- Run +1
              if(MinNRuns <= MaxNRuns){
                if(Run > MaxNRuns) {break()}
              }
              IL <- IL + PQSM2$IL[count]
              if(IL > MaxIL){
                break()}
              mRuns <- mRuns + (PQSM2$Finish[count]- PQSM2$Start[count]+1 - PQSM2$IL[count])
              if(RunFormula == TRUE){
                Formulas <- paste0(Formulas,".", PQSM2$Start[count]+1-Start,":", PQSM2$Finish[count]+1-Start)}
              if(Run >= MinNRuns & Length >= MinPQSSize){
                V1[nrow(V1)+1,] <- NA
                V1[nrow(V1),1] <- Start
                V1[nrow(V1),2] <- PQSM2$Finish[count]
                V1[nrow(V1),3] <- Run
                V1[nrow(V1),4] <- IL
                V1[nrow(V1),5] <- mRuns/Run
                if(RunFormula == TRUE){
                  V1[nrow(V1),6] <- Formulas}
              }
              count <- PQSM2$Link[count]
              if(is.na(count)== TRUE){
                break()
              }
            }}
          return(V1)
        }

        cl <- makeCluster(NCores)
        clusterExport(cl = cl, varlist = c("M2Loop","Index", "PQSM2", "RunFormula", "MaxPQSSize", "MinPQSSize", "MaxNRuns", "MinNRuns", "MaxIL"), envir = environment())
        PQSM2a <- parallel::parLapply(cl = cl, X = split(Index,seq(NROW(Index))), fun = M2Loop)
        PQSM2a <- do.call("rbind", PQSM2a)
        stopCluster(cl)

        if(RunFormula == FALSE){
          PQSM2a$Formula <- NULL
        }

      }

      PQSM2a$Length <- NULL
      PQSM2a$Sequence <- str_sub(Sequence, start = PQSM2a$Start, end = PQSM2a$Finish)
      PQSM2a$Length <- nchar(PQSM2a$Sequence)
      if(cGcC == TRUE){
        PQSM2a$cGcC <- str_sub(Sequence, ifelse(PQSM2a$Start >50, PQSM2a$Start-50, 1), ifelse(PQSM2a$Finish+50 < nchar(Sequence), PQSM2a$Finish+50 , nchar(Sequence)))
      }
      return(as.data.frame(PQSM2a))} else {PQSM2a <- data.frame()}
  }
  #
  # Detects frequency of each PQS. If Freqweight >0, it will modify final score if possible to include frequency values.
  .M2B <- function(df, RunComposition, FreqWeight){
    df$Strand <- NULL
    df$cGcC <- NULL
    df$Start <- NULL
    df$Finish <- NULL
    df <- plyr::count(df = df)
    df <- df[order(df$freq, decreasing = TRUE),]
    df <- df[,c(ncol(df), 1:(ncol(df)-1))]
    rownames(df) <- seq(1:nrow(df))
    if(!is.null(df$Score)){
      Score <- round(df$Score + ifelse(RunComposition == "C", yes = -(FreqWeight*log(df$freq,10)), no = FreqWeight*log(df$freq,10)),0)
      df$Score <- Score}
    return(df)
  }
  #
  ##### Method 3. No size MAX size limitations. No overlapping of PQS allowed.
  .M3 <- function (df, RunComposition, MinNRuns, MinPQSSize, RunFormula, Sequence, NCores, Verborrea, BulgeSize, G4hunter,  cGcC, PQSfinder){
    if(nrow(df) >= MinNRuns){
      Stops <- as.numeric(row.names(df[is.na(df$Link),]))   #Stops locations = Links that are NA :. those that are not linked
      PQSM3 <- df

      #Leader locations
      if(is.na(PQSM3$Link[1]) == TRUE){
        Leaders <- Stops[1:length(Stops)-1]+1 } else {# Leaders are the next run after a stop, except the last one
          Leaders <- c(1, Stops[1:length(Stops)-1]+1)} #seeing if the initial Run is a leader and including it as by definition it is not preceeded by a NA so it will not be detected otherwise
      if(Stops[1] == 1){Stops <- Stops[2:length(Stops)]}   #Taking down the first run as NA if it is in position 1 to force matching n? of Leaders and Stops.
      stopifnot(length(Stops) == length(Leaders))      #Number of Stops and Leaders should match.

      #Leader selectionStop
      PQSM3a <- as.data.frame(cbind(Leaders, Stops))
      Index <- PQSM3a$Leaders == PQSM3a$Stops   #Eliminating those which Leaders and stops are the same
      PQSM3a <- PQSM3a[!Index,]
      Index <- (-PQSM3a$Leaders+ PQSM3a$Stops+1 >= MinNRuns)   #Selecting those Runs that at least have MinNRuns
      PQSM3a <- PQSM3a[Index,]
      if(nrow(PQSM3a)>0){
        rownames(PQSM3a) <- seq(1:nrow(PQSM3a))

        #Building Stuff! do not disturb!
        PQSM3a$Start <- PQSM3$Start[PQSM3a$Leaders]  #Building basic parameters
        PQSM3a$Finish <- PQSM3$Finish[PQSM3a$Stops]
        PQSM3a$Length <- 1+PQSM3a$Finish -PQSM3a$Start
        PQSM3a <- PQSM3a[PQSM3a$Length >= MinPQSSize,]  #Filtering structures with at least Min desired length
        if(RunFormula == TRUE){
          PQSM3a$Formula <- paste0(PQSM3$Start[PQSM3a$Leaders]+1-PQSM3$Start[PQSM3a$Leaders],":", PQSM3$Finish[PQSM3a$Leaders]+1-PQSM3$Start[PQSM3a$Leaders])}
        rownames(PQSM3a) <- seq(1:nrow(PQSM3a))

        #Loopfunction
        M3Loop <- function(Line){
          Line <- as.numeric(Line)
          nStart <- Line[1]
          nith <- Line[2]
          Start <- Line[3]
          Finish <- Line[4]
          Length <- Line[5]

          V1 <- data.frame(Start, Finish, Length, 1, 0, 0, NA)

          count <- nStart
          V1[5] <- PQSM3$IL[nStart]
          V1[6]  <- PQSM3$Finish[nStart]-PQSM3$Start[nStart]+1- PQSM3$IL[nStart]
          if(RunFormula == TRUE){
            ForIN <- PQSM3$Start[nStart]-1
            V1[7] <- paste0(".", PQSM3$Start[nStart]-ForIN,":", PQSM3$Finish[nStart]-ForIN)
          }

          repeat {
            count <- PQSM3$Link[count]
            if(is.na(count)) {
              break()}
            V1[4] <- V1[4] + 1
            V1[5] <- V1[5] + PQSM3$IL[count]
            V1[6] <- V1[6] + (PQSM3$Finish[count]-PQSM3$Start[count]+1-PQSM3$IL[count])
            if(RunFormula == TRUE){
              V1[7] <- paste0(V1[7], ".", PQSM3$Start[count]-ForIN,":", PQSM3$Finish[count]-ForIN)}
          }
          return(V1)
        }
        if(nrow(PQSM3a)>0){
          if(NCores == 1){
            PQSM3a <- lapply(split(PQSM3a,seq(NROW(PQSM3a))),M3Loop)
            PQSM3a <- data.frame(do.call("rbind", PQSM3a))
            colnames(PQSM3a) <- c("Start", "Finish", "Length", "Runs", "IL", "mRun", "Formula")

          } else {
            cl <- makeCluster(NCores)
            clusterExport(cl = cl, varlist = c("M3Loop","PQSM3a", "PQSM3", "RunFormula"), envir = environment())
            PQSM3a <- parallel::parLapply(cl = cl, X = split(PQSM3a,seq(NROW(PQSM3a))), fun = M3Loop)
            PQSM3a <- data.frame(do.call("rbind", PQSM3a))
            colnames(PQSM3a) <- c("Start", "Finish", "Length", "Runs", "IL", "mRun", "Formula")
            stopCluster(cl)
          }

          if(RunFormula == FALSE){PQSM3a$Formula <- NULL}
          Index <- (PQSM3a$Runs >= MinNRuns)
          PQSM3a <- PQSM3a[Index,]
          PQSM3a$Sequence <- str_sub(start = PQSM3a$Start, end = PQSM3a$Finish,string =  Sequence)
          PQSM3a$mRun <- PQSM3a$mRun/ PQSM3a$Runs
          PQSM3a <- PQSM3a[order(PQSM3a$Length, decreasing = TRUE),]
          if(nrow(PQSM3a) >0){row.names(PQSM3a) <- seq(1, nrow(PQSM3a))}
          if(cGcC == TRUE){
            PQSM3a$cGcC <- str_sub(Sequence, ifelse(PQSM3a$Start >50, PQSM3a$Start-50, 1), ifelse(PQSM3a$Finish+50 < nchar(Sequence), PQSM3a$Finish+50 , nchar(Sequence)))
          }
          PQSM3a <- PQSM3a[complete.cases(PQSM3a),]
        }else {
          PQSM3a <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                             c("Start", "Finish", "Length", "Runs", "IL", "PQS", "Strand", "G4Hunter", "pqsfinder", "Score", "G", "T", "A", "C", "Conf.Quad.Seqs"))}} else {
                               PQSM3a <- setNames(data.frame(matrix(ncol = 15, nrow = 0)),
                                                  c("Start", "Finish", "Length", "Runs", "IL", "PQS", "Strand", "G4Hunter", "pqsfinder", "Score", "G", "T", "A", "C", "Conf.Quad.Seqs"))}
      return(PQSM3a)} else {PQSM3a <- data.frame()}
  }
  #
  ### Quantification of Seq
  .LoopQuantification <- function(df, LoopSeq){
    if(!is.null(df$Sequence) & length(LoopSeq) >0 & nrow(df) >0){
      df[,c(LoopSeq)] <- NA
      for(j in 1:length(LoopSeq)){
        df[,(ncol(df)-length(LoopSeq) + j)] <- round(
          stringr::str_count(df$Sequence, pattern = LoopSeq[j])*nchar(LoopSeq[j])/nchar(df$Sequence)*100, 1)
      }} else {warning(call. = T, "(The DataFrame has length 0 or hasnt the required columns) ")}
    return(df)
  }
  #
  ### Qualification
  .KnownQuadFun <- function (df, KnownNOTQuadruplex, KnownQuadruplex, RunComposition, DNA ){
    if(nrow(df) >0){
      require(dplyr)

      #Finding FUN
      KG42 <- function(RefName, RefSequence, RefLength,  df){
        Names <- rep("", nrow(df))
        Index <- NULL
        i <- which(RefLength <= df$Length)
        if(length(i)>0){
          Index <- str_detect(RefSequence, string = df$Sequence[i])
          if(sum(Index) >= 1){
            Cosa <- paste0("(?=", RefSequence, ")")
            Names[i][Index]<- paste(Names[i][Index], RefName, paste0( "(", str_count(pattern = Cosa,string =  df$Sequence[i][Index]), ")"), " ")
          }
        }
        return(Names)
      }

      if(KnownQuadruplex == TRUE){
        Ref <- filter(G4iMGrinder::Known_to_form_Quadruplex, Quadruplex == TRUE)
        if(DNA == TRUE){
          Ref <- filter(Ref, Genome == "DNA")} else {Ref <- filter(Ref, Genome == "RNA")}
        if(RunComposition == "G"){
          Ref <- filter(Ref, Nucleotide == "G")} else {Ref <- filter(Ref, Nucleotide == "C")}

        G4Name <- NULL
        AAA <- NULL

        if(nrow(Ref)>0){
          if(nrow(df) == 1){
            AAA <- mapply(FUN = KG42, RefName = Ref$Name, RefSequence= Ref$Sequence, RefLength = Ref$Length, MoreArgs = list(df[,]))
            df$Conf.Quad.Seqs  <- paste0(AAA, collapse = "")
          } else {
            if(nrow(df) > 200000){
              cof <- floor(nrow(df)/200000)
              for(i in 1:cof){
                AAA <- mapply(FUN = KG42, Ref$Name, Ref$Sequence, Ref$Length, MoreArgs = list(df[(1+((i-1)*200000)):(i*200000),]))
                G4Name <- c(G4Name, apply(AAA[,1:ncol(AAA)],1, paste0, collapse = ""))
                rm(AAA)
              }
              AAA <- mapply(FUN = KG42, Ref$Name, Ref$Sequence, Ref$Length,  MoreArgs = list(df[(1+(cof*200000)):nrow(df),]))
              G4Name <- c(G4Name, apply(AAA[,1:ncol(AAA)],1, paste0, collapse = ""))
              rm(AAA)
              df$Conf.Quad.Seqs <- G4Name
            } else {
              AAA <- mapply(FUN = KG42, RefName = Ref$Name, RefSequence= Ref$Sequence, RefLength = Ref$Length, MoreArgs = list(df[,]))
              df$Conf.Quad.Seqs  <- apply(AAA[,1:ncol(AAA)],1, paste0, collapse = "")
            }
          }
        } else {warning(call. = F, paste0("KnownQuadruplex:No data on Known-To-Form-Quadruplex fulfill those DNA/RNA and RunComposition conditions."))}
      }

      if(KnownNOTQuadruplex == TRUE){
        Ref <- filter(G4iMGrinder::Known_to_form_Quadruplex, Quadruplex == FALSE)
        if(DNA == TRUE){
          Ref <- filter(Ref, Genome == "DNA")} else {Ref <- filter(Ref, Genome == "RNA")}
        if(RunComposition == "G"){
          Ref <- filter(Ref, Nucleotide == "G")} else {Ref <- filter(Ref, Nucleotide == "C")}

        G4Name <- NULL
        AAA <- NULL

        if(nrow(Ref)>0){
          if(nrow(df) == 1){
            AAA <- mapply(FUN = KG42, RefName = Ref$Name, RefSequence= Ref$Sequence, RefLength = Ref$Length, MoreArgs = list(df[,]))
            df$Conf.NOT.Quad.Seqs  <- paste0(AAA, collapse = "")
          } else {
            if(nrow(df) > 200000){
              cof <- floor(nrow(df)/200000)
              for(i in 1:cof){
                AAA <- mapply(FUN = KG42, Ref$Name, Ref$Sequence, Ref$Length, MoreArgs = list(df[(1+((i-1)*200000)):(i*200000),]))
                G4Name <- c(G4Name, apply(AAA[,1:ncol(AAA)],1, paste0, collapse = ""))
                rm(AAA)
              }
              AAA <- mapply(FUN = KG42, Ref$Name, Ref$Sequence, Ref$Length,  MoreArgs = list(df[(1+(cof*200000)):nrow(df),]))
              G4Name <- c(G4Name, apply(AAA[,1:ncol(AAA)],1, paste0, collapse = ""))
              rm(AAA)
              df$Conf.NOT.Quad.Seqs <- G4Name
            } else {
              AAA <- mapply(FUN = KG42, RefName = Ref$Name, RefSequence= Ref$Sequence, RefLength = Ref$Length, MoreArgs = list(df[,]))
              df$Conf.NOT.Quad.Seqs  <- apply(AAA[,1:ncol(AAA)],1, paste0, collapse = "")
            }
          }
        } else {warning(call. = F, paste0("KnownNOTQuadruplex:No data on Known NOT to form quadruplex fulfill those DNA/RNA and RunComposition conditions."))}}
    }
    return(df)
  }
  #
  ### Scoring System
  # G4Hunter
  .G4Hunter <- function(df, NCores){
    if(nrow(df) > 0){
      if(NCores > 0){
        G4SCORE <- sapply(df$Sequence, .G4Hscore)
        df$G4Hunter <- round(G4SCORE*25,0)
      } else {
        #for more workers
        require(S4Vectors)
        runValue <- S4Vectors::runValue
        runLength <- S4Vectors::runLength
        Rle <- S4Vectors::Rle
        cl <- makeCluster(NCores)
        registerDoParallel(cl)
        clusterExport(cl = cl,
                      varlist = c(".G4Hscore",".G4translate", "df", "runValue", "runLength", "Rle"),
                      envir = environment())
        G4SCORE <- parSapply(cl = cl, X = df$Sequence, FUN = .G4Hscore)
        df$G4Hunter <- round(G4SCORE*25,0)
        stopCluster(cl)
      }
      #tranformed to scale of 100 to -100 from original 4 to -4
      return(df)}
  }
  # G4Hunter borrowed functions
  {

    ###### Loaned Functions for G4-Hunter from original article, Sup mat.
    #Thanks L. Lacroix
    {
      #### L. Lacroix, laurent.lacroix@inserm.fr , 20150928
      ######
      ###### G4translate change the DNA code into G4Hunter code.
      ###### Only G or C are taken into account. non G/C bases are translated in 0
      ###### It is OK if N or U are in the sequence
      ###### but might be a problem if other letters or numbers are present
      ###### lowercase ARE not welcome
      ###### G4Hscore return the G4Hscore of a sequence
      ###### The output is a number
      .G4translate <- function(x) {

        xres=x
        runValue(xres)[runValue(x)=='C' & runLength(x)>3] <- -4
        runValue(xres)[runValue(x)=='C' & runLength(x)==3] <- -3
        runValue(xres)[runValue(x)=='C' & runLength(x)==2] <- -2
        runValue(xres)[runValue(x)=='C' & runLength(x)==1] <- -1
        runValue(xres)[runValue(x)=='G' & runLength(x)>3] <- 4
        runValue(xres)[runValue(x)=='G' & runLength(x)==3] <- 3
        runValue(xres)[runValue(x)=='G' & runLength(x)==2] <- 2
        runValue(xres)[runValue(x)=='G' & runLength(x)==1] <- 1
        runValue(xres)[runValue(x)!='C' & runValue(x)!='G'] <- 0
        Rle(as.numeric(xres))
      }
      .G4Hscore <- function(y){
        # y can be DNAString or a DNAStringSet or a simple char string.

        #require(S4Vectors)
        y2 <- Rle(strsplit(as.character(y),NULL)[[1]])
        y3 <- .G4translate(y2)
        mean(y3)
      }
    }}
  # pqsfinder
  .PQSfinderfun <- function(df, RunComposition, Bt, Pb, Fm, Em, Ts, Et, Is, Ei, Ls, ET, MinNRuns){
    if(RunComposition == "G"|RunComposition == "C"){
      if(!is.null(df$Length) & !is.null(df$Runs) & !is.null(df$mRun) & nrow(df) >0){
        df$pqsfinder <- round(
          (
            (((df$mRun-1)* Bt)+Ts)^Et -                                             # Tetrad Stacking effect
              ((ifelse(is.null(df$IL), yes = 0, no = df$IL) *(Pb+1)*(MinNRuns/df$Runs))+Is)^Ei -                            # interloop in run penalizations
              ((Fm*(df$Length-ifelse(is.null(df$IL), yes = 0, no = df$IL) -(df$mRun*df$Runs))/(df$Runs-1))+Ls)^Em            # Loop penalizations
          )^Em
          ,0)
        if (RunComposition == "C") {df$pqsfinder <- df$pqsfinder*-1}
      }}
    return(df)
  }
  # cGcC
  .cGcCfun <- function(df, RunComposition){
    if(RunComposition == "G"|RunComposition=="C"){
      if(!is.null(df$cGcC) & nrow(df) >0){
        PatG <- NULL
        PatC <- NULL
        df$Gvalue <- 0
        df$Cvalue <- 0
        #foreach(icount(nrow(df)), .combine = cbind, .inorder = TRUE) %dopar% {
        for (j in 1:25){
          PatG <- paste0("(?=",paste0(rep("G", j), collapse =""),")", collapse ="")
          PatC <- paste0("(?=",paste0(rep("C", j), collapse =""),")", collapse ="")
          df$Gvalue <- df$Gvalue + (str_count(string = df$cGcC, pattern = PatG)*10*j)
          df$Cvalue <- df$Cvalue+ (str_count(string = df$cGcC, pattern = PatC)*10*j)
        }
        Index <- df$Gvalue == 0
        df$Gvalue[Index] <- 1
        Index <- df$Cvalue == 0
        df$Cvalue[Index] <- 1
        if(RunComposition != "C"){
          df$cGcC <- round(df$Gvalue/df$Cvalue,0)} else {
            df$cGcC <- round(df$Cvalue/df$Gvalue,0)*-1
          }
        df$Gvalue <- NULL
        df$Cvalue <- NULL
      }
    } else {df$cGcC <- NULL}
    return(df)}
  # Score mean
  .ScoreFun <- function(df, G4hunter, PQSfinder, cGcC, WeightParameters){
    if(nrow(df) >0){
      if (sum(WeightParameters) != 0){
        if (sum(c("G4Hunter", "pqsfinder", "cGcC") %in% colnames(df)) == 3) {
          df$Score <- round(
            (df$G4Hunter * WeightParameters[1])/sum(WeightParameters) +
              (df$pqsfinder * WeightParameters[2])/sum(WeightParameters) +
              (df$cGcC * WeightParameters[3])/sum(WeightParameters),0)
        }
        if (sum(c("G4Hunter", "pqsfinder") %in% colnames(df)) == 2 & sum(c("cGcC") %in% colnames(df))==0) {
          df$Score <- round(
            (df$G4Hunter * WeightParameters[1])/sum(WeightParameters) +
              (df$pqsfinder * WeightParameters[2])/sum(WeightParameters),0)
        }
        if (sum(c("G4Hunter", "cGcC") %in% colnames(df)) == 2 & sum(c("pqsfinder") %in% colnames(df))==0) {
          df$Score <- round(
            (df$G4Hunter * WeightParameters[1])/sum(WeightParameters) +
              (df$cGcC * WeightParameters[3])/sum(WeightParameters),0)
        }
        if (sum(c("pqsfinder", "cGcC") %in% colnames(df)) == 2 & sum(c("G4Hunter") %in% colnames(df))==0) {
          df$Score <- round(
            (df$pqsfinder * WeightParameters[2])/sum(WeightParameters) +
              (df$cGcC * WeightParameters[3])/sum(WeightParameters),0)
        }
      }}
    return(df)
  }
  # loading packages necessary
  .PackageLoading <- function(pck.CRAN = c("stringr", "stringi", "plyr", "seqinr", "stats", "parallel", "doParallel", "beepr", "stats4", "dplyr", "BiocManager"),
                            pck.BioC = c("BiocGenerics", "S4Vectors")){

    #foo was written by Simon O'Hanlon Nov 8 2013. Found in: https://stackoverflow.com/questions/4090169/elegant-way-to-check-for-missing-packages-and-install-them
    #Thanks Simon and Thanks StackOverflow and its amazing community.
    foo <- function(x){
      for( i in x ){
        #  require returns TRUE invisibly if it was able to load package
        if( ! require( i , character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)){
          #  If package was not able to be loaded then re-install
          install.packages( i , dependencies = TRUE )
          #  Load package after installing
          require( i , character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
        }
      }
    }
    foo(pck.CRAN)
    for(i in pck.BioC){
      if( ! require( i , character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE ) ){
        BiocManager::install(pck.BioC, ask = F, update = F)
      }
    }

  }
}
