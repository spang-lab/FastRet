# FastRet 1.1.0

* Added RAM caching to `getCDs

# FastRet 1.0.3

*   Initial CRAN Submission.

    Rejected because the following examples either had a "CPU time > 5s" on the CRAN testing machines or had a "CPU time > 2.5 times elapsed time". "CPU time" is the sum of the measured "user" and "system" times.

    | function             | user      | system | elapsed | ratio     |
    | -------------------- | --------- | ------ | ------- | --------- |
    | check_lm_suitability | **5.667** | 0.248  | 2.211   | **2.675** |
    | predict.frm          | 2.477     | 0.112  | 0.763   | **3.393** |
    | getCDs               | 2.745     | 0.089  | 0.961   | **2.949** |
