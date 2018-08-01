set.seed(123)

gr <- GenomicRanges::GRanges(seqnames = "Chr1",
              		     ranges = IRanges::IRanges(start = sample.int(10000, 676, replace = TRUE),
                                                       width = sample.int(10, 676, replace = TRUE),
                                                       names = paste0(rep(LETTERS[1:26], each = 26),
                                                                      rep(LETTERS[1:26], times = 26))),
                             strand = sample(c("+", "-"), size = 676, replace = TRUE))


usethis::use_data(gr, compress = "xz")

