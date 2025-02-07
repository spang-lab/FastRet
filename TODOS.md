# Todos

- [x] 1. Improve installation instructions. Done on 6.2.2025.
- [x] 2. Add Fadi as copyright holder. Done on 6.2.2025.
- [ ] 3A. Add `cdnames = CDNames` as argument to `getCDs()`, so users can decide with chemical descriptors should be returned.
- [ ] 3B. Add `cdnames = CDNames` as argument to all functions that use `getCDs()` downstream, so it can be configured from the top level calls.
- [ ] 3C. Remove all mocking options from the package and instead use `cdnames = CDNames[1:5]` in examples and tests to reduce runtime.
- [ ] 4. Make sure tests & examples don't store files outside of temp folder.
- [ ] 5. Redeploy to CRAN
- [ ] 6. Bioinformatik Teile in Paper überarbeiten
- [ ] 7. Abbildungen in Paper überarbeiten
- [ ] 8. Add visualization of feature importance
- [ ] 9. Add "Use current trained model" checkbox next to "upload trained model"
- [ ] 10. Table rows should not be selectable
- [ ] 11. Implement Fadis doc suggestions (see "..\FastRet-Hidden\Doc_Suggestions.docx")
- [ ] 12. Add myself as `contact`
- [ ] 13. Add `How to cite` section
- [ ] 14. Move `toscm/fastret-dev` to the [spang-lab container registry](https://gitlab.spang-lab.de/k8s/ha-config/container_registry)
- [ ] 15. Move `fastret-dev` and `fastret` Kubernetes deployments into HA-Config repo
- [ ] 16. Enable auto-build of `toscm/fastret-base:x.x.x` container
- [ ] 17. Show RT_ADJ ~ RT_ORIG in adjust_model plot
- [ ] 18. Enable rerender on resize
- [ ] 19. Show model coefficients
- [ ] 20. Add question mark with help texts to plots. (RS)
- [ ] 21. Move formulas from legend to help texts and use simple words in legend. (RS)
- [ ] 22. Add unit (minutes) to MSE and MAE. (RS)
- [ ] 23. Show unit in axis labels. (RS)
- [ ] 24. Rename mode "Selective Measuring" to "Design Reference Molecule Panel" (RS)
- [ ] 25. Make "minutes" as RT-unit mandatory. (RS)
- [ ] 26. Press of button "Save predictor set" currently takes multiple seconds. This should be much faster and also run asynchronous, i.e. in an "extended task".
- [ ] 27. Enable auto-build of `toscm/fastret:x.x.x-jdkyy` containers for all major jdk versions
- [ ] 28. Finalize rcdk dependency analysis
- [ ] 29. Offer "quantile normalization" as preprocessing step (RS)
- [ ] 30. Implement "remove not suitable descriptors" (Suggestion from Katja: use data from HMDB for prefiltering)
- [ ] 31. Make plots interactive (hovering over points should show molecule name)
- [ ] 32. Show PCA embedding of new points into training data after prediction
- [ ] 33. Enable auto-build of fastret-dev:x.x.x-jdkyy containers
- [ ] 34. Add a logo to the package

RS = Suggestions by Rainer Spang

# Paper

