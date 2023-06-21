source("functions.R")

library(dplyr)
library(ggplot2)
loadNamespace("GGally")
loadNamespace("caret")
loadNamespace("DMwR")
loadNamespace("ROCR")

read_data <- function(file) {
  if (tools::file_ext(file) == "feather") {
    # TODO: fallback to CSV if arrow not installed or feather file doesn't exist?
    requireNamespace("arrow", quietly = TRUE)
    df <- arrow::read_feather(file)
  } else {
    df <- read.delim(filename, header = TRUE, na.strings = c(""))
  }

  df %>%
    dplyr::filter(fluct < 2) %>% # subjects with enough data for fluctuations
    dplyr::filter(hdcat_0 == 3) %>% # subjects with manifest HD at enrollment
    select(-c("hdcat_0")) # drop constant column
}

na_counts <- function(df) {
  sapply(df, function(x) sum(is.na(x)))
}

# fluctR = 'registry', 'tfcscore'
#   'hddiagn', 'sex', 'caghigh', 'tfcscore', 'isced', 'maristat', 'res', 'fluct'
#
# fluctR-E = 'enroll', 'tfcscore'
#   'hddiagn', 'caghigh', 'tfcscore', 'region', 'isced', 'maristat', 'jobclas', 'fluct'
#
# fluctR-stage = 'registry', 'stage'
#   'hddiagn', 'sex', 'caghigh', 'stage', 'isced', 'maristat', 'res', 'fluct'
#
# fluctR-stage-E = 'enroll', 'stage'
#   'hddiagn', 'caghigh', 'stage', 'region', 'isced', 'maristat', 'jobclas', 'fluct'

list_columns <- function(name="all", stage=NULL) {
  cols <- c("hddiagn")
  if (name == "registry" | name == "all") {
    cols <- c(cols, "sex")
  }
  cols <- c(cols, "caghigh")
  if (!is.null(stage)) {
    cols <- c(cols, stage)
  }
  if (name == "enroll" | name == "all") {
    cols <- c(cols, "region")
  }
  cols <- c(cols, "isced", "maristat")
  if (name == "enroll" | name == "all") {
    cols <- c(cols, "jobclas")
  }
  if (name == "registry" | name == "all") {
    cols <- c(cols, "res")
  }
  cols <- c(cols, "fluct")
  cols
}

select_columns <- function(df, ...) {
  cols <- list_columns(...)
  df %>% dplyr::select(dplyr::all_of(cols))
}

plot_corr <- function(data) {
  GGally::ggcorr(data,
    method = c("pairwise", "spearman"),
    nbreaks = 6,
    hjust = 0.8,
    label = TRUE,
    label_size = 10,
    layout.exp = 1,
    size = 10,
    color = "grey40"
  )
}

plot_box <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_boxplot(size = 1) +
    stat_summary(
      fun = median,
      geom = "point",
      size = 3,
      color = "steelblue"
    ) +
    theme_classic(base_size = 25)
}

plot_dist <- function(data, mapping) {
  ggplot(data, mapping) +
    geom_violin() +
    theme_classic()
}

fit_model <- function(data, seed) {
  set.seed(seed)

  # tried down-sampling
  # class1<-subset(data, fluct==1)
  # class2<-subset(data, fluct==0)
  # class2<-sample_n(class2, ceiling(0.5*nrow(class1)))
  # downSample <- rbind(class2, class1)
  # rows <- sample(nrow(downSample))
  # downSample <-downSample[rows,]

  # tried up-sampling with Rose and Smote and Smote was better
  data.upscale <- DMwR::SMOTE(fluct ~ ., as.data.frame(data), perc.over = 300, perc.under = 100)

  trainDataIndex <- caret::createDataPartition(data.upscale$fluct, p = 0.75, list = F) # 75% training data
  trainData <- data.upscale[trainDataIndex, ]
  testData <- data.upscale[-trainDataIndex, ]

  # print(table(data$fluct))
  # print(table(data.upscale$fluct))
  # print(table(trainData$fluct))
  # print(table(testData$fluct))

  # use all the selected columns as ind. variables and run logit
  model <- glm(fluct ~ ., family = binomial(link = "logit"), data = trainData)
  # print(summary(model))
  # print(anova(model, test = "Chisq"))

  fitted.results <- predict(model, newdata = testData, type = "response")

  # calculate confusion matrix
  confusion <- table(testData$fluct, fitted.results > 0.5)
  # print(confusion)

  pr <- ROCR::prediction(fitted.results, testData$fluct)
  # prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  # plot(prf)

  # calculate auc
  auc <- ROCR::performance(pr, measure = "auc")
  auc <- auc@y.values[[1]]
  # print(auc)

  list(
    seed = seed,
    model = model,
    confusion = confusion,
    prediction = pr,
    auc = auc
  )
}
