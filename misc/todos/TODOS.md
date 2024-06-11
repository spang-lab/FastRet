- [x] Use async buttons for long running tasks
- [x] Make sure adjusted models are used for prediction, when available
- [x] Remove boxplots (RS)
- [x] Hide sections `Method`, `Preprocessing Options` and `Seed` behind `Show advanved Settings` button (RS)
- [x] Label "XGboost" as "XGBoost (recommended)" (RS)
- [x] Merge into main
- [x] Nach github/spang-lab umziehen
- [x] Fix missing j in "adusted"
- [x] Fix download of trained dataset
- [x] Deploy `fastret-base` to `fastret.spang-lab.de`
- [x] LATER: Move cache to `pvc`
- [x] Finalize docs incl. link to Github in pkdown site

- [ ] Set nv=3, nsw=5 in `fastret-dev`. Set maximium cores to 20 in `fastret-dev`.
- [ ] Inform Fadi

- [ ] LATER: Move `fastret-dev` to the [spang-lab container registry](https://gitlab.spang-lab.de/k8s/ha-config/container_registry)
- [ ] LATER: Enable auto-build of `toscm/fastret-base:x.x.x` container
- [ ] LATER: Show RT_ADJ ~ RT_ORIG in adjust_model plot
- [ ] LATER: Enable rerender on resize
- [ ] LATER: Show model coefficients
- [ ] LATER: Add question mark with help texts to plots. (RS)
- [ ] LATER: Move formulas from legend to help texts and use simple words in legend. (RS)
- [ ] LATER: Add unit (minutes) to MSE and MAE. (RS)
- [ ] LATER: Show unit in axis labels. (RS)
- [ ] LATER: Rename mode "Selective Measuring" to "Design Reference Molecule Panel" (RS)
- [ ] LATER: Make "minutes" as RT-unit mandatory. (RS)
- [ ] LATER: "Save predictor set" currently takes multiple seconds. This should be much faster and also run asynchronous, i.e. in an "extended task".
- [ ] LATER: Enable auto-build of `toscm/fastret:x.x.x-jdkyy` containers for all major jdk versions
- [ ] LATER: Finalize rcdk dependency analysis
- [ ] LATER: Offer "quantile normalization" as preprocessing step (RS)
- [ ] LATER: Implement RAM caching for `getCD`
- [ ] LATER: Implement "remove not suitable descriptors" (Suggestion from Katja: use data from HMDB for prefiltering)
- [ ] LATER: Make plots interactive (hovering over points should show molecule name)
- [ ] LATER: Show PCA embedding of new points into training data after prediction
- [ ] LATER: Enable auto-build of fastret-dev:x.x.x-jdkyy containers

RS = Suggestions by Rainer Spang