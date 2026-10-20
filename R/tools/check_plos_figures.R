# Check figure files against the PLOS requirements for TIFF figures, the
# checks that PLOS's PACE tool makes, before uploading them.
#
# Usage: Rscript R/tools/check_plos_figures.R <folder or .tif files ...>
#        e.g. the TIFFs of stage 80 (results/figures/plos/) and Fig 5 from the
#        screenshot tool of the app repository
#
# For each file, reads the TIFF header and checks:
#   - one image per file (no further pages or layers)
#   - LZW compression
#   - RGB or greyscale, 8 bits per channel, no alpha channel
#   - 300 to 600 dpi
#   - width 2.63 to 7.5 in, height at most 8.75 in
#   - file size at most 10 MB
# Prints one line per file (OK, or the problems found) and exits with status 1
# if any file fails. Nothing is written.

args <- commandArgs(TRUE)
if (length(args) == 0) stop("usage: Rscript R/tools/check_plos_figures.R <folder or .tif files ...>")
files <- unlist(lapply(args, function(a) {
  if (dir.exists(a)) list.files(a, pattern = "\\.tiff?$", full.names = TRUE, ignore.case = TRUE) else a
}))
if (length(files) == 0) stop("no .tif files found in: ", paste(args, collapse = " "))

# --- TIFF header ------------------------------------------------------------

# the tags of the first image (IFD), and the number of images in the file
read_tiff_tags <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  bytes <- readBin(con, "raw", file.size(path))
  order <- rawToChar(bytes[1:2])
  if (!order %in% c("II", "MM")) stop("not a TIFF file")
  endian <- if (order == "II") "little" else "big"
  uint <- function(offset, size) {
    readBin(bytes[(offset + 1):(offset + size)], "integer", size = size, signed = size == 4, endian = endian)
  }
  if (uint(2, 2) != 42) stop("not a classic TIFF file")
  type_size <- c(1, 1, 2, 4, 8)  # BYTE, ASCII, SHORT, LONG, RATIONAL
  value <- function(type, count, field) {
    size <- type_size[type] * count
    at <- if (size <= 4) field else uint(field, 4)
    switch(as.character(type),
      "3" = vapply(seq_len(count) - 1, function(i) uint(at + 2 * i, 2), 1),
      "4" = vapply(seq_len(count) - 1, function(i) uint(at + 4 * i, 4), 1),
      "5" = uint(at, 4) / uint(at + 4, 4),
      NA)
  }
  ifd <- uint(4, 4)
  tags <- list()
  n_images <- 0
  while (ifd != 0) {
    n_images <- n_images + 1
    n <- uint(ifd, 2)
    if (n_images == 1) {
      for (i in seq_len(n) - 1) {
        entry <- ifd + 2 + 12 * i
        tags[[as.character(uint(entry, 2))]] <- value(uint(entry + 2, 2), uint(entry + 4, 4), entry + 8)
      }
    }
    ifd <- uint(ifd + 2 + 12 * n, 4)
  }
  list(tags = tags, n_images = n_images)
}

# --- checks -----------------------------------------------------------------

check_figure <- function(path) {
  tiff <- tryCatch(read_tiff_tags(path), error = function(e) e)
  if (inherits(tiff, "error")) return(conditionMessage(tiff))
  tag <- function(code, default = NA) {
    x <- tiff$tags[[as.character(code)]]
    if (is.null(x)) default else x
  }
  width <- tag(256); height <- tag(257)
  compression <- tag(259, 1); photometric <- tag(262)
  samples <- tag(277, 1); bits <- tag(258, 1)
  dpi <- tag(282); unit <- tag(296, 2)
  if (!is.na(dpi) && unit == 3) dpi <- dpi * 2.54  # resolution per cm
  problems <- c(
    if (tiff$n_images > 1) sprintf("%d images or layers, expected 1", tiff$n_images),
    if (compression != 5) sprintf("compression %d, expected LZW (5)", compression),
    if (!photometric %in% c(1, 2)) sprintf("colour model %d, expected RGB (2) or greyscale (1)", photometric),
    if (!is.null(tiff$tags[["338"]]) || samples > if (identical(photometric, 2)) 3 else 1) "alpha or extra channel",
    if (any(bits != 8)) sprintf("%s bits per channel, expected 8", paste(unique(bits), collapse = "/")),
    if (is.na(dpi) || unit == 1) "no resolution (dpi) recorded",
    if (!is.na(dpi) && (dpi < 299.5 || dpi > 600.5)) sprintf("%.0f dpi, expected 300 to 600", dpi),
    if (!is.na(dpi) && width / dpi > 7.5 + 1e-6) sprintf("%.2f in wide, at most 7.5", width / dpi),
    if (!is.na(dpi) && width / dpi < 2.63) sprintf("%.2f in wide, at least 2.63", width / dpi),
    if (!is.na(dpi) && height / dpi > 8.75 + 1e-6) sprintf("%.2f in high, at most 8.75", height / dpi),
    if (file.size(path) > 10e6) sprintf("%.1f MB, at most 10", file.size(path) / 1e6)
  )
  size <- if (is.na(dpi)) sprintf("%d x %d px", width, height) else
    sprintf("%.2f x %.2f in, %.0f dpi, %.1f MB", width / dpi, height / dpi, dpi, file.size(path) / 1e6)
  if (length(problems) == 0) paste("OK", size) else paste(c(size, problems), collapse = "; ")
}

results <- vapply(files, check_figure, "")
for (f in files) cat(sprintf("%-14s %s\n", basename(f), results[[f]]))
failed <- !startsWith(results, "OK")
if (any(failed)) {
  cat(sprintf("\n%d of %d files do not meet the PLOS requirements.\n", sum(failed), length(files)))
  quit(status = 1)
}
cat(sprintf("\nAll %d files meet the PLOS requirements.\n", length(files)))
