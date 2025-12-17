# Package Internals

## Model Training

Model training is implemented in the `train_model` function. It first
tries out several different parameter values in cross validation,
remembers the best parameter value and then uses this value to train a
final model on the full dataset. For further details see the function
documentation of
[train_model](https://spang-lab.github.io/FastRet/reference/train_frm.html).

## Reactive Graph

FastRet is a reactive web app, that means most of its calculations are
triggered as reaction to user inputs. The reactive graphs for FastRet’s
four modes are shown below, with

- Action buttons shown as red rectangles
- Upload buttons shown as orange rectangles
- Other inputs shown as green rectangles
- Intermediate calculations shown as violet rectangles
- Outputs shown as blue rectangles
- Synchronous functions shown as white rectangles with two lines on the
  left and right side
- Asynchronous functions shown as hexagons
- Red arrows from one node to another indicating that the first node
  triggers the second node
- Green arrows from one node to another indicating that the first node
  holds an input value for the second node
- Blue arrows from one node to another indicating that the first node
  writes an output value to the second node

### Train new Model

[![train_new_model.svg](Package-Internals/train_new_model.svg)](https://spang-lab.github.io/FastRet/articles/Package-Internals/train_new_model.svg)

### Predict Retention Times

[![predict_retention_times.svg](Package-Internals/predict_retention_times.svg)](https://spang-lab.github.io/FastRet/articles/Package-Internals/predict_retention_times.svg)

### Selective Measuring

[![selective_measuring.svg](Package-Internals/selective_measuring.svg)](https://spang-lab.github.io/FastRet/articles/Package-Internals/selective_measuring.svg)

### Adjust Existing Model

[![adjust_existing_model.svg](Package-Internals/adjust_existing_model.svg)](https://spang-lab.github.io/FastRet/articles/Package-Internals/adjust_existing_model.svg)
