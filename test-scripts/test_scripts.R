devtools::load_all()


## Test sharing index function -------------------------------
peak_files <- list.files("../../Projects/insulinomas/data/H3K27ac/peaks",
                         pattern="^NET.*\\.broadPeak",
                         full.names=TRUE)

bam_files <- list.files("../../Projects/insulinomas/data/H3K27ac/bams",
                        pattern="^NET.*\\.bam$",
                        full.names=TRUE)

ids <- c("NET17", "NET29")

peak_files <- peak_files[sapply(ids, grep, peak_files)]
bam_files <- bam_files[sapply(ids, grep, bam_files)]

plot_sharing_index(peak_files=peak_files,
                   bam_files=bam_files,
                   plot_title="Insulinoma H3K27ac")

## Test conservation scores ------------------------------
peak_file <- "../../Projects/insulinomas/data/H3K27ac/peaks/HI_19_peaks.broadPeak"

# peak_file <- regioneR::toGRanges(peak_file)[1:1e4,]
# colnames(mcols(peak_file))[1] <- "test_changing_id"

cons <- get_conservation_scores(peak_file)

ggplot(cons, aes(pos, phastCons)) +
  # geom_ribbon(aes(group=type, ymin=phastCons-sd, ymax=phastCons+sd)) +
  geom_line(aes(color=type), lwd=1)
