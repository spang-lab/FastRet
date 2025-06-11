test_that("RP_Mod$Normal should be a strict subset of RP", {
  mod_idx <- match(RP_Mod[[1]]$SMILES, RP$SMILES)
  expect_false(any(is.na(mod_idx)))
  RP_Sub <- RP[mod_idx, c("NAME", "SMILES", "RT")]
  rownames(RP_Sub) <- NULL
  expect_identical(RP_Sub, RP_Mod$Normal)
  # Because RP_Mod$Normal only contains metabolites that are present in RP and
  # we ensure that all datasets in RP_Mod have the same metabolites (next test),
  # we also know that all other datasets in RP_Mod are subsets of RP.
})

test_that("All datasets in RP_Mod should have the same metabolites", {
  RP_Mod_1 <- RP_Mod[[1]][, c("NAME", "SMILES")]
  for (i in names(RP_Mod)) {
    RP_Mod_i <- RP_Mod[[i]][, c("NAME", "SMILES")]
    expect_identical(RP_Mod_i, RP_Mod_1)
  }
})

test_that("No metabolite of RP_Val[[i]] should be included in RP", {
  expect_false(any(RP_Val[[1]]$Name %in% RP$Name))
  expect_false(any(RP_Val[[1]]$SMILES %in% RP_Mod[[1]]$SMILES))
  # It's sufficient to test this for the first dataset in RP_Val, because in the
  # next test we ensure that all datasets in RP_Val have the same metabolites.
})

test_that("All datasets in RP_Val should have the same metabolites", {
  RP_Val_1 <- RP_Val[[1]][, c("NAME", "SMILES")]
  for (i in names(RP_Val)) {
    RP_Val_i <- RP_Val[[i]][, c("NAME", "SMILES")]
    expect_identical(RP_Val_i, RP_Val_1)
  }
})
