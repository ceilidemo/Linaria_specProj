## 01_PLSDA : exploratory PLSDA on first few time pts
# Ceili DeMarais

library(dplyr)
library(caret)
library(klaR)
library(agricolae)
library(corrplot)
library(reshape)
library(ggplot2)

#setup
datasets <- list(
  raw = list(dat = "data_work/spec_clean.csv", work = "data_work/PLSDA/raw"),
  processed = list(dat = "data_work/spec_clean_processed.csv", work = "data_work/PLSDA/processed")
)

for (dataset_name in names(datasets)) {
  
  dat  <- datasets[[dataset_name]]$dat
  work <- datasets[[dataset_name]]$work
  dir.create(work, recursive = TRUE, showWarnings = FALSE)
  
  cat("\n\n", dataset_name, "\n")
  
  dati <- read.csv(dat) %>%
    mutate(
      treatment = factor(treatment),
      leafPos   = factor(leafPos),
      timePt    = factor(timePt)
    ) %>%
    dplyr::filter(leafPos == "bot") %>%
    mutate(group = paste(treatment, leafPos, sep = "_"))
  
  wvl_cols    <- names(dati)[grepl("^X[0-9]", names(dati))]
  wvl_seq     <- wvl_cols[seq(1, length(wvl_cols), by = 10)]
  time_points <- levels(dati$timePt)
  
  ###################################################
  ## PART 1 : Quick PLSDA to determine kappas #######
  ###################################################
  for (tp in time_points) {
    
    cat("\nSample:", tp, "\n")
    
    dati_tp <- dati %>% dplyr::filter(timePt == tp)
    classi  <- factor(dati_tp$group)
    spec    <- dati_tp[, wvl_seq]
    
    n_per_class <- min(table(dati_tp$group))
    n_train     <- round(n_per_class * 0.70)
    cat("n per class:", n_per_class, "| n_train:", n_train, "\n")
    
    nsims <- 10
    rndid <- list()
    set.seed(1840)
    for (i in 1:nsims) {
      rndid[[i]] <- with(dati_tp, ave(1:nrow(dati_tp),
                                      group,
                                      FUN = function(x) {sample.int(length(x))}))
    }
    
    mods <- list()
    for (nsim in seq(nsims)) {
      inTrain    <- rndid[[nsim]] <= n_train
      print(nsim)
      flush.console()
      set.seed(nsim)
      traini     <- spec[inTrain, ]
      testi      <- spec[!inTrain, ]
      trainclass <- classi[inTrain]; testclass <- classi[!inTrain]
      plsFit <- train(traini, trainclass, method = "simpls", tuneLength = 20,
                      trControl = trainControl(method = "boot"))
      mods[[nsim]] <- plsFit
    }
    
    ncomps <- vector(length = nsims)
    for (i in 1:nsims) ncomps[i] <- mods[[i]]$finalModel$ncomp
    table(ncomps)
    
    kappas_tune <- data.frame(ncomps = 1:20, matrix(NA, nrow = 20, ncol = length(mods)))
    for (i in 1:length(mods)) kappas_tune[, i + 1] <- mods[[i]]$results$Kappa
    
    kapp <- as.data.frame(as.numeric(t(kappas_tune[, -1])))
    kapp <- cbind(kapp, rep(1:20, each = length(mods)))
    names(kapp) <- c("Kappa", "ncomps")
    kapp$ncomps <- as.factor(kapp$ncomps)
    
    modi    <- lm(Kappa ~ ncomps, kapp)
    tuk     <- HSD.test(modi, "ncomps")
    tuk_dat <- as.data.frame(tuk$groups)
    tuk_dat$var <- as.numeric(row.names(tuk_dat))
    tuk_dat <- tuk_dat[order(tuk_dat$var, decreasing = F), ]
    letters <- as.character(tuk_dat$groups)
    
    pdf(file.path(work, paste0(tp, "_kappas.pdf")), width = 5, height = 4)
    par(bty = "l")
    boxplot(kapp$Kappa ~ kapp$ncomps, ylim = c(0, 1),
            xlab = "Number of components", ylab = "Kappa",
            main = paste(tp, "-", dataset_name))
    text(x = 1:20, y = rep(1, 20), letters)
    dev.off()
    
  } 
} 
  
  
###################################################
## PART 2 : FULL RUN! #############################
###################################################

#determine comps per dataset based on their kappa plot
compi_per_dataset <- list(
  raw = list(
    data_0 = 7,
    data_1 = 7,
    data_2 = 2,
    data_3 = 7,
    data_4 = 3
  ),
  processed = list(
    data_0 = 5,
    data_1 = 5,
    data_2 = 5,
    data_3 = 7,
    data_4 = 3
  )
)

for (dataset_name in names(datasets)) {
  
  dat  <- datasets[[dataset_name]]$dat
  work <- datasets[[dataset_name]]$work
  
  cat("\n\n", dataset_name, "\n")
  
  dati <- read.csv(dat) %>%
    mutate(
      treatment = factor(treatment),
      leafPos   = factor(leafPos),
      timePt    = factor(timePt)
    ) %>%
    dplyr::filter(leafPos == "bot") %>%
    mutate(group = paste(treatment, leafPos, sep = "_"))
  
  wvl_cols    <- names(dati)[grepl("^X[0-9]", names(dati))]
  wvl_seq     <- wvl_cols[seq(1, length(wvl_cols), by = 10)] #take every 10th wvl
  time_points <- levels(dati$timePt)
  
  for (tp in time_points) {
    
    cat("\nSample:", tp, "\n")
    
    compi <- compi_per_dataset[[dataset_name]][[tp]]
    cat("compi:", compi, "\n")
    
    dati_tp <- dati %>% dplyr::filter(timePt == tp)
    classi  <- factor(dati_tp$group)
    spec    <- dati_tp[, wvl_seq]
    
    n_per_class <- min(table(dati_tp$group))
    n_train     <- round(n_per_class * 0.70)
    
    nsims <- 50 #change this based on #iterations you want
    rndid <- list()
    set.seed(1840)
    for (i in 1:nsims) {
      rndid[[i]] <- with(dati_tp, ave(1:nrow(dati_tp),
                                      group,
                                      FUN = function(x) {sample.int(length(x))}))
    }
    
    finmods <- list()
    for (nsim in 1:nsims) {
      print(nsim)
      flush.console()
      set.seed(nsim)
      inTrain    <- rndid[[nsim]] <= n_train
      training   <- spec[inTrain, ]
      testing    <- spec[!inTrain, ]
      trainclass <- as.factor(classi[inTrain]); testclass <- as.factor(classi[!inTrain])
      
      finalModel <- plsda(training, trainclass, ncomp = compi, probMethod = "Bayes", method = "simpls")
      finmods[[nsim]] <- finalModel
    }
    
    probis <- list()
    confus <- list()
    
    for (nsim in seq(nsims)) {
      print(nsim)
      flush.console()
      set.seed(nsim)
      inTrain    <- rndid[[nsim]] <= n_train
      testing    <- spec[!inTrain, ]
      testclass  <- as.factor(classi[!inTrain])
      
      plsProbs   <- predict(finmods[[nsim]], newdata = testing, type = "prob")
      plsClasses <- predict(finmods[[nsim]], newdata = testing)
      confus[[nsim]] <- confusionMatrix(data = plsClasses, testclass)
      
      probs <- as.data.frame(plsProbs)
      names(probs) <- sapply(strsplit(names(probs), split = ".n"), "[", 1)
      probs <- cbind(testclass, probs)
      probis[[nsim]] <- probs
    }
    
    # Prob plot (didnt actually make one... i dont like them)
    arr       <- array(unlist(probis), dim = c(dim(probis[[1]]), nsims))
    prob_mean <- apply(arr, 1:2, mean)
    prob_mean <- as.data.frame(prob_mean)
    prob_mean$V1 <- probis[[1]]$testclass
    colnames(prob_mean) <- colnames(probis[[1]])
    
    pp <- melt(prob_mean, id = "testclass")
    pp$position  <- ifelse(as.character(pp$testclass) == as.character(pp$variable), 2, 1)
    pp$testclass <- factor(pp$testclass, levels = rev(levels(pp$testclass)))
    
    # Confusion matrix
    tabs  <- list()
    for (i in 1:length(confus)) tabs[[i]] <- confus[[i]]$table
    tabsi    <- Reduce("+", tabs)
    tab_mean <- as.data.frame.matrix(tabsi / length(confus))
    sums     <- colSums(tab_mean)
    
    tabs_perc <- matrix(NA, length(sums), length(sums))
    for (i in 1:length(sums)) tabs_perc[, i] <- tab_mean[, i] / sums[i]
    colnames(tabs_perc) <- colnames(confus[[1]]$table)
    rownames(tabs_perc) <- rownames(confus[[1]]$table)
    
    col <- colorRampPalette(c("black", "black", "brown", "gold", "forestgreen"))
    
    pdf(file.path(work, paste0(tp, "_confumatrix.pdf")), width = 7, height = 6, pointsize = 11)
    corrplot(tabs_perc, p.mat = tabs_perc, insig = "p-value", sig.level = -1, addCoef.col = 1,
             tl.srt = 70, col = col(20), cl.lim = c(0, 1), tl.col = 1, tl.offset = 1.5,
             cl.ratio = 0.2, cl.align.text = "l", cl.cex = 0.9, mar = c(1, 3, 3, 3))
    mtext("Prediction", 2, line = 0, cex = 1.2)
    mtext("Reference", at = 2, line = 1.5, cex = 1.2)
    dev.off()
    
    # Loadings
    lls <- list()
    for (i in 1:length(finmods)) {
      lls[[i]] <- abs(loadings(finmods[[i]])[1:dim(loadings(finmods[[1]]))[1], 1:compi])
      sumis <- lapply(lls, rowSums)
    }
    
    mm <- apply(simplify2array(sumis), 1, mean)
    ss <- apply(simplify2array(sumis), 1, sd)
    mm <- as.data.frame(mm)
    mm <- cbind(mm, ss)
    
    pdf(file.path(work, paste0(tp, "_loadings.pdf")), width = 6, height = 4)
    plot(1:ncol(testing), mm$mm, type = "n", bty = "l", ylab = "abs(loadings)",
         xaxt = "n", xlab = "Wavelength (nm)", main = paste(tp, "-", dataset_name),
         ylim = c(min(mm$mm) - mm[which(mm$mm == min(mm$mm)), ]$ss,
                  max(mm$mm) + mm[which(mm$mm == max(mm$mm)), ]$ss))
    polygon(x = c(1:ncol(testing), ncol(testing):1),
            y = c(mm$mm - mm$ss, rev(mm$mm + mm$ss)), col = "grey", border = "grey")
    lines(1:ncol(testing), mm$mm, type = "l")
    axis(1, at = seq(1, nrow(mm), length.out = 6), labels = seq(400, 2400, 400))
    dev.off()
    
    # Stats
    specifi <- list(); sensi <- list(); preci <- list()
    mean_specifi <- list(); mean_sensi <- list()
    accu   <- numeric(length = length(confus))
    kappas <- numeric(length = length(confus))
    
    for (i in 1:nsims) {
      confu   <- confus[[i]]$table
      n       <- sum(confu)
      nc      <- nrow(confu)
      diag_v  <- diag(as.matrix(confu))
      rowsums <- apply(confu, 1, sum)
      colsums <- apply(confu, 2, sum)
      p <- rowsums / n
      q <- colsums / n
      preci[[i]]  <- diag_v / colsums
      accuracy    <- sum(diag_v) / n
      accu[i]     <- accuracy
      expAccuracy <- sum(p * q)
      kappas[i]   <- (accuracy - expAccuracy) / (1 - expAccuracy)
      sensi[[i]]  <- diag_v / rowsums
      speci <- numeric(length = nrow(confu))
      for (j in 1:nrow(confu)) speci[j] <- sum(confu[-j, -j]) / sum(confu[, -j])
      names(speci)      <- colnames(confu)
      specifi[[i]]      <- speci
      mean_specifi[[i]] <- sum(speci) / nc
      mean_sensi[[i]]   <- sum(sensi[[i]]) / nc
    }
    
    modelstatistics_table <- data.frame(
      Metric = c("Sensitivity", "Specificity", "Precision", "Accuracy", "Kappa"),
      Mean   = c(mean(unlist(sensi)), mean(unlist(specifi)), mean(unlist(preci)),
                 mean(accu), mean(kappas)),
      SD     = c(sd(unlist(sensi)), sd(unlist(specifi)), sd(unlist(preci)),
                 sd(accu), sd(kappas))
    )
    
    write.csv(modelstatistics_table, file.path(work, paste0(tp, "_modelstatistics.csv")), row.names = FALSE)
    cat("Accuracy:", round(mean(accu), 3), "| Kappa:", round(mean(kappas), 3), "\n")
    
  }
}
