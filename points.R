library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(Lithics3D)
library(dplyr)

setwd("./Tabun")

# Create list of files with .pts extension
filelist <- list.files(pattern = ".pts")

# Get names from files
names <- gsub(".pts", "", filelist)

# Create empty array
landmarks <- NULL

# Import each file and paste them in a single matrix
for (i in 1:length(filelist)){
  tmp <- as.data.frame(read.table(filelist[i], skip = 2, header = F)[,2:4])
  fix <- tmp[1:3,]
  tmp <- arrange(tmp, desc(row_number()))
  tmp <- rbind(fix,tmp)
  tmp <- tmp[1:63,]
  tmp <- as.matrix(tmp)
  landmarks <- rbind(landmarks, tmp)
}

landmarks <- mapply(landmarks, FUN=as.numeric)
landmarks <- matrix(data=landmarks, ncol=3)

# Convert landmark 3D data matrix into array
landmarks <- arrayspecs(landmarks, 63, 3)

# Spread landmarks on all specimens for GM (Will's approach)
lms <- array(0, dim=c(40,3,49))

for(i in 1:dim(landmarks)[3]) {
  pp <- landmarks[1,,3]
  startA <- landmarks[3,,3]
  startB <- landmarks[2,,3]
  A <- digit.curves(startB, rbind(landmarks[59:63,,3], pp), nPoints = 3, closed = F)
  B <- digit.curves(pp, rbind(landmarks[44:48,,3], startA), nPoints = 3, closed = F)
  C <- digit.curves(startA, rbind(landmarks[49:58,,3]), nPoints = 7, closed = F)
  D <- digit.curves(startA,rbind(landmarks[4:43,,3],startB), nPoints = 30, closed = F)
  comb <- rbind(pp, startB, startA, A[2:(nrow(A)-1),], B[2:(nrow(B)-1),], C[2:(nrow(C)-1),], D[2:(nrow(D)-1),])
  lms[,,i] <- comb
}


# Same thing but to use in Elliptic Fourier Analysis (difference is that a single curve is computed)
for(i in 1:dim(landmarks)[3]) {
  pp <- landmarks[1,,i]
  startA <- landmarks[3,,i]
  startB <- landmarks[2,,i]
  slms_all <- digit.curves(startA,rbind(landmarks[4:43,,i],startB,landmarks[59:63,,i],pp,landmarks[44:48,,i], startA), nPoints = 28, closed = F)
  lms[,,i] <- slms_all
}


# Plot lmks

specimen1 <- vcgPlyRead("Point1933.790.B.1_LR.ply")
specimen1 <- specimen1 <- vcgUniformRemesh(specimen1, voxelSize=median(vcgMeshres(specimen1)$edgelength))

shade3d(specimen1, color = "gray")
spheres3d(lms[1:2,,1], color = "red", size = 10)
spheres3d(lms[3:32,,1], color = "blue", size = 10)
spheres3d(lms[33:40,,1], color = "green", size = 10)



# Procrustes

points_gpagen <- geomorph::gpagen(lms)
PCA <- gm.prcomp(points_gpagen$coords)
plot(PCA)

msh <- mshape(points_gpagen$coords)
plotRefToTarget(PCA$shapes$shapes.comp1$min, msh)
plotRefToTarget(PCA$shapes$shapes.comp1$min, PCA$shapes$shapes.comp1$max, method= "vector", mag = )

scores <- as.data.frame(PCA$x)

plot(scores[,1],points_gpagen$Csize)



landmarks[1,,] <- x %>% arrange(desc(row_number()))

###############
### Create atlas



atlas <- createAtlas(specimen1, landmarks = comb[1:46,], patch = pat)
plotAtlas(atlas)


patched <- placePatch(atlas = atlas, dat.array = specimen2, path = ".", inflate = 5, prefix = c("specimen2"), fileext = ".ply")


checkLM(patched, path=".", atlas=atlas, suffix = "")



