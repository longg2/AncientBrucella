library(fastbaps)
library(ape)

sparse.data <- import_fasta_sparse_nt("ST11FullSeq.fasta")
sparse.data <- optimise_prior(sparse.data, type = "optimise.symmetric")
multi <- multi_res_baps(sparse.data)
