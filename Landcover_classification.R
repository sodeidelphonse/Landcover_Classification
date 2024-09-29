
#==================
# GROUP 2 PROJECT
#==================


#==============================================================================
# PART 1: PREDICTION OF LAND COVER CLASSES
#==============================================================================

library(terra)

# 1. load the croped sentinel data from yesterday
sentinel_may = rast("data/S2_20230501_ndvi.tif")
sentinel_nov = rast("data/S2_20221127_ndvi.tif")

# Combine sentinel scenes
sentinel =  c(sentinel_may, sentinel_nov)
names(sentinel)
plot(sentinel)

# load the polygons: can also be a .shp file
reference_data <- vect("data/Reference_data.gpkg")
head(reference_data, 3)
table(reference_data$class)

plotRGB(sentinel,
  b = "B02_may",
  g = "B03_may",
  r = "B04_may",
  stretch = "lin"
)
plot(reference_data, add = TRUE, col = "red")


# 2. Extract spectral data for training polygons ===============================

df <-  extract(sentinel, reference_data, na.rm = TRUE)
reference_data$ID <- seq(1, nrow(reference_data))
df <- merge(df, reference_data)


# 3. Train a Random Forest model ==============================================
library(randomForest)

rfmodel = randomForest(x = df[,c("B02_may", "B03_may", "B04_may", "B08_may", "NDVI_may",
                                 "B02_nov", "B03_nov", "B04_nov", "B08_nov", "NDVI_nov")],
                      y = as.factor(df[, c("class")]),
                      ntree = 100)


# 4. Apply the model to the image -- classification ===========================

# Apply the model to the sentinel image
lcc <-  predict(sentinel, rfmodel)
plot(lcc)

plot(lcc, col = c("darkgreen", "grey10", "lightblue", "grey40", "firebrick1", "lightgreen", "limegreen"))

writeRaster(lcc, "data/outputs/landcover_model1.tif", overwrite = TRUE)

land_cover_classes <- data.frame(levels(lcc))
write.csv(land_cover_classes, "data/outputs/land_cover_classes.csv", row.names = FALSE)


## Cross validation ===========================================================

# 3. Split reference data into training (70% training data, 30% test data)
library(caret)

set.seed(125)

trainIndex <- createDataPartition(df$class, p = 0.7, list = FALSE)

# Split data into training and testing 
trainData <- df[trainIndex, ]
head(trainData, 5)

testData  <- df[-trainIndex, ]
head(randomForest)

# 4. Apply the Random Forest model
rfmodel <- randomForest(x = trainData[,c("B02_may", "B03_may", "B04_may", "B08_may", "NDVI_may", "B02_nov", "B03_nov", "B04_nov", "B08_nov", "NDVI_nov")],
                       y = as.factor(trainData[, c("class")]),
                       ntree = 500)

# Apply the model to the sentinel image
lcc <- predict(sentinel, rfmodel)
plot(lcc, col = c("darkgreen", "grey10", "lightblue", "grey40", "firebrick1", "lightgreen"))

# 5. Perform an accuracy assessment of the independent test set
test_prediction <-  predict(rfmodel, testData)

# Confusion matrix
cfm <- table(testData$class, test_prediction)
cfm

# accuracy: sum of diagonals divided by all
sum(diag(cfm)/sum(cfm))


# 6. use k-fold validation  ====================================================
## Cross-validation

#set.seed(124)
# set up cross validation with 5 folds 
trc = trainControl(method = "cv", number = 5)

# tune model mtry
rfmodel_cv_may = caret::train(x = trainData[,c("B02_may", "B03_may", "B04_may", "B08_may", "NDVI_may")],
                              y = trainData[, c("class")],
                              method = "rf",
                              trControl = trc,
                              ntree = 100,
                              tuneLength = 3)

rfmodel_cv_may

rfmodel_cv_nov = caret::train(x = trainData[,c("B02_nov", "B03_nov", "B04_nov", "B08_nov", "NDVI_nov")],
                              y = trainData[, c("class")],
                              method = "rf",
                              trControl = trc,
                              ntree = 100,
                              tuneLength = 3)

rfmodel_cv_nov

rfmodel_cv_full = caret::train(x = trainData[,c("B02_may", "B03_may", "B04_may", "B08_may", "NDVI_may", "B02_nov", "B03_nov", "B04_nov", "B08_nov", "NDVI_nov")],
                               y = trainData[, c("class")],
                               method = "rf",
                               trControl = trc,
                               ntree = 100,
                               tuneLength = 3)

lcc2 <- predict(sentinel, rfmodel_cv_full)
plot(lcc2, col = c("darkgreen", "grey10", "lightblue", "grey50", "firebrick1", "lightgreen"))

# Predict the classes of the independent test set
test_prediction <- predict(rfmodel_cv_full, testData)

# Confusion matrix
cfm <- confusionMatrix(as.factor(testData$class), test_prediction, mode = "everything")

write.csv(cfm$table, "data/outputs/confusion_matrix.csv", row.names = FALSE)
write.csv(cfm$byClass, "data/outputs/accuracy_by_class.csv", row.names = FALSE)
write.csv(cfm$overall, "data/outputs/accuracy_overall.csv", row.names = FALSE)

# Prediction = user's accuracy
# Recall = Producer's accuracy

writeRaster(lcc2, "data/outputs/landcover_model2.tif", overwrite = TRUE)

land_cover_classes <- data.frame(levels(lcc2))
write.csv(land_cover_classes, "data/outputs/land_cover_classes2.csv", row.names = FALSE)


#==============================================================================
# Part 2: Are there any difference in terms of amount of IS among land covers
#==============================================================================

library(dplyr)
library(ggplot2)

#... Load Fogo field data

fogo_species <- read.csv("data/fogo_species.csv")
head(fogo_species)
str(fogo_species)

fogo_plots <- read.csv("data/fogo_plots.csv")
head(fogo_plots)
str(fogo_plots)


#... Data processing and analysis =============================================

# join the field data
fogo_merged <- merge(fogo_plots, fogo_species, by = "PLOT")
head(fogo_merged)
str(fogo_merged)

# the introduced species diversity and amount
prop_spec <- length(table(fogo_merged$SPECIES[fogo_merged$SOURCE =="IS"]))/length(table(fogo_merged$SPECIES))
prop_spec  # 57.34% of the observed single species

prop_amount <- sum(table(fogo_merged$SPECIES[fogo_merged$SOURCE =="IS"]))/sum(table(fogo_merged$SPECIES))
prop_amount 

# convert the field data into a spatial vector object
plot_points <- vect(fogo_merged, 
                    geom = c("LONG", "LAT"), 
                    crs = "EPSG:4326")  

# project the data into UTM system
crs(sentinel)
plot_points <- project(plot_points, "EPSG:32626")

# Extract land cover classes at species locations
df_lcc      <- extract(lcc2, plot_points)
df_classes <- cbind(df_lcc, plot_points, geom(plot_points)[, c("x", "y")])
head(df_classes)

# Save the results
write.csv(df_classes, "data/combined_data.csv")


#... Distribution of species amount ==========================================

# Extract the introduced species (IS) information 
spec_df <- df_classes %>% 
  filter(SOURCE == "IS") %>%
  select(PLOT, x, y, class, SPECIES) %>%
  group_by(PLOT, x, y, class) %>%
  summarise(amount = n())
head(spec_df)

# Distribution of introduced species amount
pal1 <- c("darkgreen", "grey10", "grey50", "firebrick1", "lightgreen")
ggplot(spec_df, aes(x = class, y = amount, fill = class)) +
  geom_boxplot() +
  scale_fill_manual(values = pal1, guide = "none") +
  #theme_minimal() +
  labs(x = "Landcover", y = "Species amount") 

# The most eight abundant introduced species
df_classes %>% 
  filter(SOURCE == "IS") %>%
  group_by(SPECIES) %>% 
  summarise(amount = n()) %>% 
  arrange(desc(amount)) %>%
  head(8)

# Test the significant difference of IS amounts among land cover classes
head(spec_df)
kruskal.test(spec_df$amount ~ spec_df$class)
# P < 0.001 (there is a significant difference in terms of IS amount)


#... Distribution of IS on the map ============================================

#cap <- vect("cabo_verde/gadm41_CPV_1.shp")
#cap <- project(cap, "EPSG:32626")

# plot of land cover and sample sites
pal2 <- c("darkgreen", "grey10", "lightblue", "grey50", "firebrick1", "lightgreen")

plot(lcc2, col = pal2)
#plot(cap, add = TRUE)
points(plot_points[plot_points$SOURCE == "IS",], col = "white", cex = 0.5) #, add = TRUE)

## NB: a complete map can be designed using GIS tools or ggplot2/tmap packages


#==============================================================================
##  END OF THE SCRIPT
#==============================================================================

