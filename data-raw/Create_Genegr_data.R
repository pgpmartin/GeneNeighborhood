set.seed(123)

Genegr <- GenomicRanges::GRanges(seqnames = "Chr1",
               		         ranges = IRanges::IRanges(start = sample.int(10000, 676, replace = TRUE),
                                                           width = sample.int(10, 676, replace = TRUE),
                                                           names = paste0(rep(LETTERS[1:26], each = 26),
                                                                          rep(LETTERS[1:26], times = 26))),
                                 strand = sample(c("+", "-"), size = 676, replace = TRUE))

GenomeInfoDb::seqinfo(Genegr) <- GenomeInfoDb::Seqinfo(seqnames = "Chr1",
                                                       seqlengths = 10010,
                                                       isCircular = FALSE,
                                                       genome = "mock")

usethis::use_data(Genegr, compress = "xz")

