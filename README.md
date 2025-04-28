# plier2

## Comparing PLIERv2 with PLIERv1 (original)

1. Create a conda environment for PLIERv2, and _another_, different conda
   environment for PLIERv1.
1. Run PLIERv2 following `gtex.Rmd`. The input GTEx file can be downloaded from
   [here](https://olucdenver-my.sharepoint.com/:f:/r/personal/marc_subiranagranes_cuanschutz_edu/Documents/plier2?csf=1&web=1&e=gAAyQl).

   Once you reach the line at the end that runs `PLIERv2`, wrap those lines with this so we can measure the time it takes:
   ```r
   start.time <- Sys.time()

   # Run PLIERv2
   gtex.plier=PLIERv2(gtexFBMfiltered, matchedPaths, sdres = gtex.sdres, Chat = Chat,doCrossval = T, multiplier = 3, max.U.updates = 3, max.iter = 100)

   end.time <- Sys.time()
   time.taken <- end.time - start.time
   time.taken
   ```
1. Once PLIERv2 finishes running, save the PLIERv2 parameteres to external `rds` files so you can load them in
   the PLIERv1 environment:
   ```r
   dir.create("example", showWarnings = FALSE)
   saveRDS(gtexFBMfiltered[], "example/data.rds")
   saveRDS(as.matrix(matchedPaths), "example/allpaths.rds")
   saveRDS(gtexFBM.svd, "example/svdres.rds")
   saveRDS(as.matrix(Chat), "example/chat.rds")
   ```
1. Load those files into your PLIERv1 session and run PLIER:
   ```r
   data<-readRDS("example/data.rds")
   allpaths<-readRDS("example/allpaths.rds")
   svdres<-readRDS("example/svdres.rds")
   chat<-readRDS("example/chat.rds")

   library(PLIER)

   start.time <- Sys.time()

   # Run PLIERv1
   plierResult=PLIER(data, allpaths, k=200, svdres=svdres, Chat=chat, max.iter=100)

   end.time <- Sys.time() 
   time.taken <- end.time - start.time
   time.taken
   ```
