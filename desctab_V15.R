
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#% Name und Speicherort des Programms (Pfad): desctab.R
#% Projekt/Studie: keine
#% Ersterstellung (Autor, Datum): Michael Kramer, 17.07.2015
#% Zweck des Programms: Tabelle mit deskriptiven Statistiken erstellen
#% Modifikationen (Autor, Datum, Beschreibung): welcher Lageparameter und welches
#%     Streuungsmaß bei stetigen Variablen ausgegeben wird, wird jetzt aus dem 
#%     Rowname ausgelesen (z.B. Alter, median (range): median und range werden ausgelesen,
#%     ermittelt und ausgegeben)
#%     Michael Kramer, 05.11.2015, Medianes Survival mit IQR/CI eingebaut
#%     Michael Kramer, 21.12.2015, Tabellierung der kategorialen Variablen so umgebaut, dass
#%        immer Spalte "TRUE" ausgegeben wird und missings immer mitgezählt werden
#%     Michael Kramer, 03.02.2016, Header fett ausgeben, eigene sanitize.text.function, damit 
#%        subtitle Zeilen fett formatiert werden können
#%     Michael Kramer, 02.03.2016, option useNA hinzugefügt, um bei kategorialen Variablen NAs
#%        variable behandeln zu können; Ermittlung des 25%-Quartil bei IQR des Survival korrigiert
#%     Michael Kramer, 11.03.2016, Flexibilisierung der Ergebnisdarstellung durch Verwendung von gsub
#%     Michael Kramer, 15.03.2016, added p-value adjustment
#%     Michael Kramer, 22.08.2016, Fehler korrigiert, der bei der Bildung der Einflussvariable für die
#%        Signifikanztests bei der Definition von subsets dafür sorgte, dass auch irrelevante Gruppen mit
#%        getestet wurden; fehlende Werte bei Chisq-Test standardmäßig ausschließen
#%     Michael Kramer, 25.08.2016, separates Argument useNA.test="no" eingeführt um missings
#%        in der Tabellendarstellung und im Test separat behandeln zu können, Standardwert für useNA wieder auf
#%        "ifany" gesetzt, um standardmäßig die Prozente der Gesamtzahl angezeigt zu bekommen, und nicht nur die der
#%        nichtfehlenden Werte
#%     Michael Kramer, 26.08.2016, grundlegende Überarbeitung: Zerlegung in Funktion für Schätzer und Funktion 
#%        für Tests; Steuerung, was zu tun ist über Taskmatrix
#%     Michael Kramer, 26.08.2016, subset-Option zur Analyse der Häufigkeiten von Subsets und 
#%        korrekte Anzeige von Prozenten für die Subsets eingefügt, subtitle in varset umbenannt
#%     Michael Kramer, 31.08.2016, suboff keyword eingeführt, um den Bereich für Subset-Analysen besser 
#%        begrenzen zu können
#% Version: V15-0
#% Vorlage: V14-0
#% Input (z.B. SAS-Dateien, Makrovariablen): 
#%  depvar, gelabelte Liste mit abhängigen Variablen, Listenname erscheint in der Spalte "Parameter" der 
#%    ausgegenen Tabelle. Beispiel:
#%    depvar <- list("Age (years), median (range)"=d$AufAlter,"Male gender, n (%)"=d$Sex=="m  ",
#%                   "subtitle"="Blutbild"); Angabe der gewünschten Lage- und Streuungsparameter
#%    nach einem Komma definiert, welche Parameter ausgegeben werden (mean=Mittelwert;
#%    median=Median; SD=Standardabweichung; IQR=Interquartilbereich; range=Minimum-Maximum); 
#%    normal range=2,5% Quantil bis 97,5% Quantil
#%    mittels "subtitle" kann eine Zwischenüberschrift in die Tabelle eingefügt werden, im Beispiel steht
#%    dann "Blutbild" in der Spalte "Parameter" und die restlichen Spalten bleiben leer
#%  invar, gelabelte Liste mit Gruppen, für die Statistiken für depvar erstellt werden sollen
#%    Gruppen müssen als "logical" oder "factor" angegeben sein
#%    invar <- list("All patients n = 2200"=1,"AML96 n = 1862"=d$TRIALID=="AML96")
#%  digits, gibt an auf wiviele Nachkommastellen Ergebnisse gerundet werden sollen (Integer-Wert angeben)
#% Output (z.B. SAS-Dateien, Makrovariablen): data.frame mit den Spalten Parameter und names(invar) und Zeilen
#%    names(depvar)
#% Namen von Programmen, 
#% welche vorher gelaufen sein m?ssen: keine
#% Namen von Programmen, 
#% welche im Programm aufgerufen werden: keine
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


desctab <- function(invar=NULL,depvar,test=FALSE,digits=2,ndigits=0,
                    useNA="ifany",useNA.test="no",
                    padjust=NULL){
  # Formatierungen für Latex
  options(xtable.sanitize.text.function = function(x){gsub("%","\\\\%",x)})
  
  # Classes of dependent variables
   dclass <- unlist(lapply(depvar,class))            # Variablenklassen in taskmatrix eintragen
  
  ###########################
  # Für subsetting Länge der 
  # Variablenvektoren in depvars
  ###########################
  # Auswahl der Listenelemente, die Daten enthalten
   sel <- which(!grepl("multinom",names(depvar))&
                !grepl("subset",names(depvar))&
                !grepl("subtitle",names(depvar)))
                      
  # Variablenlängen bestimmen
   varlen <- unlist(lapply(depvar[sel],length))
  # Survival-Objekte sind doppelt so lang, deswegen muss die varlen hier halbiert werden
   varlen[which(dclass[sel]=="Surv")] <- varlen[which(dclass[sel]=="Surv")]/2
  # unique varlen generieren und überprüfen, ob sie genau 1 ist
   varlen <- unlist(unique(varlen))
  
  # warnings and errors
   if(length(varlen)>1) stop("lengths of depvars differ")
   
  ###########################
  # if an element of invar has lenght 1
  # statistics should be calcualted for
  # all observations, to do so a logical
  # vector full of 'TRUE' with varlen
  # is generated and substitutes the 
  # invar-element with length 1
  ###########################
   ilen <- lapply(invar,length)
   if(1%in%ilen) invar[[which(ilen==1)]] <- rep(TRUE,varlen)
  
  # warnings and errors
   invarlen <- unlist(unique(lapply(invar,length)))
   if(length(invarlen)>1) stop("lengths of invars differ")
   if(!invarlen==varlen) stop("lengths of invar and depvar differ")
   
  ###########################  
  # Task-Matrix generieren
  ###########################
  # Taskmatrix extrahiert die Argumente aus den Listennamen und ordnet
  # sie in einer Matrix an, auf die dann bei der Asuführung zugrgriffen
  # wird, damit das Programm weiß, was jetzt zu tun ist

  # extract names from depvar
    nam <- names(depvar) # Labels der labelled list
  
  # find and indicate keywords
    multinom <- grepl("multinom",nam)     # start of variable varid for multinomial categorical variable
    multioff <- grepl("multioff",nam)     # end of variable varid for multinomial categorical variable
    varid <- rep(NA,length(depvar))
    subset <- grepl("subset",nam)         # Definition of start of subset
    suboff <- grepl("suboff",nam)         # Definition end of subset
    subsetid <- rep(NA,length(depvar))    # Welches Listenelement bildet das subset?
    subtitle <- grepl("subtitle",nam)     # ist das Element ein subtitle

  # remove keywords
    nam <- gsub("subset","",nam)
    nam <- gsub("suboff","",nam)
    nam <- gsub("multinom","",nam)
    nam <- gsub("multioff","",nam)
    nam <- gsub("subtitle","",nam)
  
  # write subsetid into all rows belonging to the same subset
  # rows of a subset are all rows from row with keyword 'subset'
  # to row with keywort 'suboff'
    substart=FALSE
    for(i in 1:length(subset)){
      if( subset[i])                         {subsetid[i] <- i            ; substart=TRUE}
      if(!subset[i] & !suboff[i] & substart)  subsetid[i] <- subsetid[i-1]
      if(!subset[i] &  suboff[i] & substart) {subsetid[i] <- subsetid[i-1]; substart=FALSE}
    }
    subsetid[which(subset)] <- NA
  
  # write varid into all rows belonging to the same multinomial variable
    multistart=FALSE
    for(i in 1:length(multinom)){
      if( multinom[i])                             {varid[i] <- i         ; multistart=TRUE}
      if(!multinom[i] & !multioff[i] & multistart)  varid[i] <- varid[i-1]
      if(!multinom[i] & multioff[i] & multistart)  {varid[i] <- varid[i-1]; multistart=FALSE}
      if(!multinom[i] & !multioff[i] & !multistart) varid[i] <- i
    }
    nvarid <- rep(table(varid),table(varid))          # besteht eine Gruppe aus mehr als einer Zeile?

  # warnings and errors
   if(!table(multinom)["FALSE"]==table(multioff)["FALSE"]) stop("multinom started, but not ended")
   if(!table(subset)["FALSE"]==table(suboff)["FALSE"]) stop("subset started, but not ended")
  
  
  # write all estimates in all rows where estimates should be printed in the
  # results table
    # extract estimates from nam
     charnum <- unlist(gregexpr(",",nam))+1            # Estimate vorhanden?
     estimates <- substr(nam,charnum+1,nchar(nam))     # Estimates extrahieren
     estimates[charnum==0] <- NA
     es <- !is.na(estimates)                 # Indicator for definition of an estimate
    # write estimates in rows from estimate defining rows 
    # into where estimates are not given
     if(length(depvar)>1){
       for(i in 2:length(estimates)){
          if(is.na(estimates[i])&!subset[i]) estimates[i] <- estimates[i-1] 
       }
     }
    estimates[multinom==TRUE|subtitle==TRUE] <- NA
    estimates <- gsub("\\}","",estimates)  # remove geschweifte brackets
  

  
  # indicate rows, where test results should be printed in the results table
  # a test should be done for each variable were an estimate is desired, or
  # one for a group of variables indicated by 'multinom'
  # no test should be calculated for categorical variables, where no percentage is
  # desired
    if(test){
     testlocal <- es
     testlocal[!is.na(estimates) & nvarid==1] <- TRUE
     testlocal[estimates=="n"] <- FALSE
    }else{
      testlocal <- rep(FALSE,length(depvar))
    }
     
  # combine the vectors to the taskmatrix
    tm <- data.frame(nam,subtitle,multinom,multioff,varid,subset,suboff,subsetid,dclass,estimates,testlocal) 
  
  # nam       = name of parameter and estimate, when given
  # multinom  = indicator for definition of a multinomial variable in this row, estimate of this row
  #             is used for all consequtive rows, until keyword 'multioff' is recognized, row 
  #             with 'multioff' is included as last row in the multinomial-variable set
  # multioff  = indicates the end of a multinomial variable; row with the keyword 'multioff' is
  #             the last element in the multinomial variable
  # varid     = indicates which rows build one multinomial variable in the depvar-list
  # subset    = indicator for definition of subset of data, the statistics will be applied to;
  #             an estimate can be desired, but doesn't have to; subset definition is valid until
  #             keyword 'suboff' is recognized in one of the following rows; row with keyword 'suboff'
  #             is the last row subset is valid for
  # suboff    = indicates the end for subset analysis; row with keyword 'suboff' is the last element
  #             that will be analysed in the subset
  # subsetid  = rownumber of depvars, where the subset is defined for the next analyses,
  #             if is.na all observations will be used
  # estimates = estimates that is desired for the row, if is.na no estimate will be calcualted
  # dclass    = class of dependent variable in the row
  # testlocal = if test==TRUE, indicates, that a test is calculated and printed in the row where
  #             testlocal == TRUE
  
  #############################
  # initialize results data.frame
  #############################
   out <- matrix(NA,nrow=length(depvar),ncol=length(names(invar))) 
   out <- data.frame(out)
   # Spalten der Ergebnismatrix benennen
   names(out) <- c(names(invar))

  
  
  #############################
  # invar is given as separate TRUE/FALSE
  # indicator variables; to use standard
  # test functions of R a single variable
  # with as much levels as groups are
  # defined is built here; because it is
  # possible to define groups in invar,
  # that are not distint it is tested 
  # whether groups are distinct. If they
  # are not distinct, no tests will be 
  # calculated, because argument test is 
  # set to FALSE
  #############################
  
  if(test==TRUE){
    # Gruppenvariable generieren, weil diese in die separaten invars aufgeteilt ist
     testivar <- droplevels(interaction(rbind.data.frame(invar)))
     testivar[grepl("TRUE",testivar)==FALSE] <- NA
     testivar <- droplevels(testivar)
    # Wenn Gruppen nicht distinct, dann test auf FALSE zurücksetzen
     tiv <- rbind.data.frame(invar)
     tval <- max(as.numeric(unlist(apply(tiv,1,sum))),na.rm=TRUE) # tval > 1, wenn Gruppen nicht distinct
     test <- tval==1
    if(length(invar)<2) test <- FALSE # nur eine Gruppe
    if(test==FALSE)  warning("Groups are not distinct, no tests will be calculated")
    if(test==TRUE) out$"p.value" <- NA
  }
  
  #############################
  # Funktion für Test
  #############################
  
  do.test <- function(depvars,testivar,...){
    # Test durchführen und p-Wert ausgeben
    # Kategoriale Variable
    if(tm$dclass[i]%in%c("logical","factor")) {
      # Testvariable generieren, weil diese in separate depvars aufgeteilt sein kann
      if(length(depvars)>1){    # varset Zeile ohne Daten entfernen
        depvars <- depvars[-1]
      }
           
      if(length(depvars)>1){    # wenn mehrere Variablen, dann per interaction zu einer kombinieren
        testdvar <- droplevels(interaction(rbind.data.frame(depvars)))
        testdvar[grepl("TRUE",testdvar)==FALSE] <- NA
        testdvar <- droplevels(testdvar)
      }else{
        testdvar <- unlist(depvars)
      }
      # Wenn Gruppen nicht distinct, dann test auf FALSE zurücksetzen
      tdv <- rbind.data.frame(depvars)
      tval <- max(as.numeric(unlist(apply(tdv,1,sum))),na.rm=TRUE) # tval > 1, wenn Gruppen nicht distinct
      test <- tval==1
      
      if(test==TRUE){
        tab <- table(testivar,testdvar,useNA=useNA.test)
        p <- try(round(chisq.test(tab)$p.value,4))
        # Fehlermeldung bei fehlgeschlagenem Chisq-Test entfernen
        if(!is.numeric(p)) p <- NA
      }else{
        p <- NA
      }
    }
    # Stetige Variable
    if(tm$dclass[i]%in%c("numeric","integer")) {
      p <- try(round(kruskal.test(unlist(depvars)~testivar)$p.value,4))
      # Fehlermeldung bei fehlgeschlagenem KW-Test entfernen
      if(!is.numeric(p)) p <- NA
    }
    # Survival Variable
    if(tm$dclass[i]=="Surv"){
      dvar <- depvars[[1]]
      p <- try(round(pchisq(survdiff(dvar~testivar)$chisq,
                            df=length(levels(testivar))-1,
                            lower.tail=FALSE),4))
      # Fehlermeldung bei fehlgeschlagenem KW-Test entfernen
      if(!is.numeric(p)) p <- NA
    }
    return(p)
  }
  
  #############################
  # Funktion für Schätzer
  #############################
  
  do.est <- function(dvar,ivar,...){
    est <- tm$estimates[i]
    dvar <- unlist(dvar)
    ivar <- unlist(ivar)
    if(length(ivar)==1) ivar <- rep(TRUE,length(dvar))
    # numerische Abhängige Variable
    if(tm$dclass[i]%in%c("numeric","integer")){
      # Schätzer ermitteln
      smry <- summary(dvar[ivar==TRUE])
      est <- gsub("mean",round(smry["Mean"],digits),est)
      est <- gsub("median",round(smry["Median"],digits),est)
      est <- gsub("SD",round(sd(dvar[ivar==TRUE],na.rm=TRUE),digits),est)
      # Streuungsmaß ermitteln
      est <- gsub("range",paste(round(smry["Min."],digits),
                                round(smry["Max."],digits),sep="-"),est)
      est <- gsub("IQR",paste(round(smry["1st Qu."],digits),
                              round(smry["3rd Qu."],digits),sep="-"),est)
    }
    
    # kategoriale Abhängige Variable
    if(tm$dclass[i]%in%c("logical","factor")){
      # absolute Häufigkeit ermitteln
      tab <- table(dvar,ivar,useNA=useNA)
      # keine leeren Zellen in Kreuztabelle
      if("TRUE"%in%rownames(tab) & 
           "TRUE"%in%colnames(tab)){
        est <- gsub("n",tab["TRUE","TRUE"],est) 
        # Prozente ermitteln
        ptab <- prop.table(tab,2)
        est <- gsub("%",round(ptab["TRUE","TRUE"]*100,ndigits),est) 
      }
      # leere Zellen in Kreuztabelle
      if(!"TRUE"%in%rownames(tab)|!"TRUE"%in%colnames(tab)){
        est <- gsub("n",0,est)
        est <- gsub("%",0,est)
      }
    }
    
    # Survival Variable als abhängige
    if(tm$dclass[i]=="Surv"){
      if(!is.na(table(ivar)["TRUE"])){ # gültige Fälle vorhanden
        # Schätzer ermitteln
         kmfit <- survfit(dvar[ivar==TRUE]~1)
         est <- gsub("median",round(min(kmfit$time[kmfit$surv<0.5]),digits),est)
        # IQR ermitteln
         est <- gsub("IQR",paste(round(min(kmfit$time[kmfit$surv<0.75]),digits),
                            round(min(kmfit$time[kmfit$surv<0.25]),digits),sep="-"),est)
        # Konfidenzintervall ermitteln
         est <- gsub("95%-CI",paste(round(min(kmfit$time[kmfit$lower<0.5],na.rm=TRUE),digits),
                               round(min(kmfit$time[kmfit$upper<0.5],na.rm=TRUE),digits),sep="-"),est)
         est <- gsub("CI",paste(round(min(kmfit$time[kmfit$lower<0.5],na.rm=TRUE),digits),
                           round(min(kmfit$time[kmfit$upper<0.5],na.rm=TRUE),digits),sep="-"),est)
      }else{
        est <- "NA"
      }
    }
    return(est)
  }
  
  
  
  
  #############################
  # Variablen abarbeiten
  #############################
  for(i in 1:length(depvar)){ # loop um abhängige Variable
    # build subset, if defined via subset keyword
    if(!is.na(tm$subsetid[i])){
      sg <- depvar[[tm$subsetid[i]]] # für subset, alle mit TRUE
      if(!TRUE%in%sg) stop(paste("empty subset defined in line",i-1))
    }
    # no subset defined
    if(is.na(tm$subsetid[i])){
       sg <- rep(TRUE,varlen) # für alle
    }
    
    for(j in 1:length(invar)){  # loop um Einflussvariable
      # calculate test
       if(tm$testlocal[i]==TRUE & j==1){ # only calculate test in 1st run for i
         out$"p.value"[i] <- do.test(depvars=lapply(depvar[which(varid==tm$varid[i])], "[", which(sg)),testivar[which(sg)])
       }
      # calculate statistics
       if(!is.na(tm$estimates[i])){
         out[i,j] <- do.est(dvar=depvar[[i]][sg],ivar=invar[[j]][sg])
       }
    }
  }
  out <- cbind("Parameter"=nam,out)
  
  if(!is.null(padjust)){
    out$padj <- p.adjust(out$"p.value",method=padjust)
  }
  return(out)
}