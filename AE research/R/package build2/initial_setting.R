#working_dir <- "/home/proj/jeffwu/isye-ae"
working_dir <- "Y:/AE test2"
design_file_name <- "design.csv"
rawdata_file_path <- "Y:/AE test2/RawReactingData"

Case.Index <- 1:5
Response.ColumnIndex <- 7      # select flow responses
Response.ColumnName <- c("T")
Snapshot.Index <- 0:199          # select the snapshop
SnapshotForPOD.Index <- 0:199     # select the snapshop for POD
if(any(!is.element(SnapshotForPOD.Index, Snapshot.Index)))
  stop("SnapshopForPOD.Index Setting is Wrong!!")
Rescale_by_1000.ColumnIndex <- 6  # select the response index that needs to be rescaled
split.num <- 5                    # split the inner product matrix: 
                                  # if the inner product matrix is too big, suggest using this setting
POD_Mode_Energy.cutoff <- 0.99    # POD mode energy
POD_Mode_number.cutoff <- 200     # the most number of POD modes