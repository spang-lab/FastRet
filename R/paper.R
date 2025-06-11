reproduce_paper <- function() {
    # VALIDATE: number of metabolites in in-house RT library (~380)
    # VALIDATE: model adjustment works with as few as 25 metabolites
    # VALIDATE: BRT MAE for RP (0.5), HILIC (0.86), RP-AX (X.XX)
    # VALIDATE: min-size for metabolites for adjustment
    # VALIDATE: number of metabolites in each dataset (Retip HILIC: 970, In-house HILIC: 364, Mixed mode: 412, RP: 40)
    # VALIDATE: RMSE, MAE, R2, PPB1M for all datasets and methods
    # VALIDATE: 200 candidate lambda values, 10-fold CV for Lasso
    # VALIDATE: xgboost hyperparameters (tree depth 2-5, nrounds 300-1000 step 100, eta 0.01/0.02), 10-fold CV
    # VALIDATE: selective measuring method as described is implemented and used (steps 1-5)
    # VALIDATE: k = 25 is used for medoid selection
    # VALIDATE: 25 medoids are used for further validation
    # VALIDATE: adjustment is performed using polynomial regression (linear, quadratic, cubic) on medoids
    # VALIDATE: 10-fold cross-validation is used for all performance evaluations
    # VALIDATE: BRTs outperform Lasso by 0.26 RMSE and 0.24 MAE on average
    # VALIDATE: PPB1M for RP (85.95%), RP-AX (46.37%)
    # VALIDATE: HILIC in-house (364 metabolites, 327 per fold), Retip (970, 873 per fold), in-house HILIC performs slightly better
    # VALIDATE: Figure 1 shows RMSE vs training set size for RP, saturation point
    # VALIDATE: Lasso saturates at ~25, BRTs at ~60 training samples
    # VALIDATE: For small sets Lasso > BRTs, for large sets BRTs > Lasso
    # VALIDATE: external validation set contains 25 metabolites not in training
    # VALIDATE: Figure 2 shows RT changes for all tested conditions
    # VALIDATE: original model MAE = 1 min on validation set
    # VALIDATE: adjusted model MAE = 0.65 min
    # VALIDATE: medoids-only model MAE = 1.06 min, original model trained on 401 metabolites
    # VALIDATE: adjusted models outperform medoids-only and original model
    # VALIDATE: Table values match actual MAE for each approach
    # VALIDATE: performance comparable to Retip paper
    # VALIDATE: in-house HILIC slightly better than Retip HILIC
    # VALIDATE: BRTs > Lasso for large datasets
    # VALIDATE: Lasso reliable for small training sets
    # VALIDATE: chemical descriptor quality is limiting factor (if possible to check)
}
