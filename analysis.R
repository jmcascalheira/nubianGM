library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(Lithics3D)

setwd("./3D/")

# Create list of files with .pts extension
filelist <- list.files(pattern = ".pts")

# Get names from files
names <- gsub(".pts", "", filelist)

# Create empty array
landmarks <- NULL

# Import each file and paste them in a single matrix
for (i in 1:length(filelist)){
  tmp <- as.matrix(read.table(filelist[i], skip = 2, header = F)[,2:4])
  tmp <- tmp[71:120,]
  landmarks <- rbind(landmarks, tmp)
}


# Convert landmark 3D data matrix into array
landmarks <- arrayspecs(landmarks, 50, 3)

# Sliding all specimens

lms <- array(0, dim=c(30,3,49))

for(i in 1:dim(landmarks)[3]) {
  startA <- landmarks[1,,i]
  startB <- landmarks[40,,i]
  slms_A <- digit.curves(startA, landmarks[2:40,,i], nPoints = 20, closed = F)
  slms_B <- digit.curves(startB, rbind(landmarks[40:50,,i],startA), nPoints = 8, closed = F)
  comb <- rbind(startA, startB, slms_A[2:(nrow(slms_A)-1),], slms_B[2:(nrow(slms_B)-1),])
  lms[,,i] <- comb
}





# Plot lmks

specimen1 <- vcgPlyRead("...")

shade3d(specimen1)
points3d(lms[1:2,,1], color = "red", size = 10)
points3d(lms[3:22,,1], color = "blue", size = 10)
points3d(lms[23:30,,1], color = "green", size = 10)



# Procrustes

twoD <- lms[,1:2,]


flake_gpagen <- gpagen(lms)
PCA <- gm.prcomp(flake_gpagen$coords)
plot(PCA)

msh <- mshape(flake_gpagen$coords)
plotRefToTarget(PCA$shapes$shapes.comp1$min, msh)
plotRefToTarget(PCA$shapes$shapes.comp1$min, PCA$shapes$shapes.comp1$max, method = "vector", mag = 2)











arrange(desc(row_number()))


