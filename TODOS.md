# Planned

## Categorize existing Excel sheets

List all files in `misc/datasets` and add a description to each one.

## Add "Use current trained model" checkbox

Add "Use current trained model" checkbox next to "upload trained model"

## Add a logo to the package

## Fill out NEWS.md entries before 1.0.3

Required to list contributions in paper.

## Improve Examples

- Add `cdnames = CDNames` as argument to `getCDs()`, so users can decide
    witch chemical descriptors should be returned.
- Add `cdnames = CDNames` as argument to all functions that use
    `getCDs()` downstream, so it can be configured from the top level calls.
- Remove all mocking options from the package and instead use `cdnames =
    CDNames[1:5]` in examples and tests to reduce runtime.

## Move DockerHub images to thespanglab

Steps:

1. Push an improved example image to Dockerhub/thespanglab
2. Delete the fastret-dev deployment from Kubernetes
3. Assert the fastret deployment uses an image from gitlab.spang-lab.de
4. Delete the FastRet image from Dockerhub/toscm
5. Check that the example image is mentioned in the FastRet documentation

Optional:

1. Move the FastRet image to the Containers group
2. Update the FastRet image in the FastRet deployment

Super optional:

1. Move the Lirec image to the Containers group
2. Update the Lirec image in the Lirec deployment

## Bioinformatik Teile in Paper überarbeiten

## Abbildungen in Paper überarbeiten

## Table rows should not be selectable

## Add myself as contact

## Add How to cite section

## Move fastret-dev image to HA-config

- Move `toscm/fastret-dev` image from DockerHub to the [spang-lab container registry](https://gitlab.spang-lab.de/k8s/ha-config/container_registry)
- Move `fastret-dev` and `fastret` Kubernetes deployments into HA-Config repo

## Add visualization of feature importance

## Implement Fadis doc suggestions

See "..\FastRet-Hidden\Doc_Suggestions.docx"

## Show RT_ADJ vs. RT_ORIG in adjust_model plot

## Enable rerender on resize

## Show model coefficients

## Add help texts to plots

Add question mark with help texts to plots. Suggested by Rainer.

## Move formulas from legend to help text

Move formulas from legend to help texts and use simple words in legend. Suggested by Rainer.

## Show that MSE and MAE is given in minutes

## Show unit in axis labels

Suggested by Rainer.

## Make minutes as RT-unit mandatory

Suggested by Rainer.

## Speedup "Save predictor set"

Press of button "Save predictor set" currently takes multiple seconds.
This should be much faster and also run asynchronous, i.e. in an "extended task".

## Enable auto-build of FastRet containers

- Enable auto-build of `toscm/fastret-base:x.x.x` container
- Enable auto-build of `toscm/fastret:x.x.x-jdkyy` containers for all major jdk versions

## Make plots interactive (hovering over points should show molecule name)

## Show PCA embedding of new points into training data after prediction


# Done

## Improve installation instructions

Done on 6.2.2025.

## Add Fadi as copyright holder

Done on 6.2.2025

## Redeploy to CRAN

Done around March 2025.

## Always use temp in tests and examples

Done around March 2025.
