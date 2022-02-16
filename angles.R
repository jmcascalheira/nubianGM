library(Morpho)
library(Rvcg)
library(geomorph)
library(rgl)
library(Lithics3D)


# Set file location (only if needed)
setwd("./3D")

# Import model
specimen1 <- vcgPlyRead("XXXX.ply")

# Remesh model to make more uniform for the calculation of angles (this might take some time to run)
specimen1 <- vcgUniformRemesh(specimen1, voxelSize=median(vcgMeshres(specimen1)$edgelength))


# Read landmarks
landmarks <- as.matrix(read.table("CoreME78.627Aa_LR_distal_LMs.pts", skip = 2, header = F)[,2:4])


# Get path using landmarks and following ridges
edge.curve.ids <- sPathConnect(landmarks[84:86,], specimen1, path.choice="ridges")
edge.curve.path = t(specimen1$vb)[edge.curve.ids,1:3]


# Plot specimen, landmarks and line to check
shade3d(specimen1, color = "green")
points3d(landmarks[84:86,], color="purple", size=10)
lines3d(edge.curve.path, color="blue", size=6, lwd=4)

# Distribute 5 equally spaced points across the line
edge.curve = pathResample(edge.curve.path, 5, method="npts")

# Export table with 5 landmarks
write.table(x = edge.curve[2:6,], file = "equal_space_coord.xyz", row.names = F, col.names = F)

# Calculate angles for the 3 mid landmarks (extremes are not calculated due to the possible interference of other surfaces)
res = edgeAngles(specimen1, edge.curve, 3)


library(dplyr)
angles <- dplyr::select(data.frame(t(sapply(res,c))), angle)
angles <- apply(angles,2,as.character)
write.table(x = angles, file = "angles.txt", row.names = F, col.names = F)



rgl.close()

# Plot landmarks and lines for each angle 
points3d(edge.curve, color="purple", size=10)

for(i in 1:length(res)){
  lines3d(res[[i]]$inters.pts, lwd=3)
  ang.l = rbind(res[[i]]$inters.pts[1,],
                edge.curve[i+1,],
                res[[i]]$inters.pts[2,])
  lines3d(ang.l, color="blue", lwd=1)
}



