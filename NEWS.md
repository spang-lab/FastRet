# FastRet 1.1.0

* Added RAM caching to `getCDs()`

# FastRet 1.0.3

*   Initial CRAN Submission.

    Rejected because the following examples caused at least one of the following notes on the CRAN testing machines: `CPU time > 5s`, `CPU time > 2.5 times elapsed time`. In this context, `CPU time` is calculated as the sum of the measured `user` and `system` times.

    | function             | user      | system | elapsed | ratio     |
    | -------------------- | --------- | ------ | ------- | --------- |
    | check_lm_suitability | **5.667** | 0.248  | 2.211   | **2.675** |
    | predict.frm          | 2.477     | 0.112  | 0.763   | **3.393** |
    | getCDs               | 2.745     | 0.089  | 0.961   | **2.949** |
