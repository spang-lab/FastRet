# Start the FastRet GUI

Starts the FastRet GUI

## Usage

``` r
start_gui(port = 8080, host = "0.0.0.0", reload = FALSE, nw = 2, nsw = 1)
```

## Arguments

- port:

  The port the application should listen on

- host:

  The address the application should listen on

- reload:

  Whether to reload the application when the source code changes

- nw:

  The number of worker processes started. The first worker always
  listens for user input from the GUI. The other workers are used for
  handling long running tasks like model fitting or clustering. If `nw`
  is 1, the same process is used for both tasks, which means that the
  GUI will become unresponsive during long running tasks.

- nsw:

  The number of subworkers each worker is allowed to start. The higher
  this number, the faster individual tasks like model fitting can be
  processed. A value of 1 means that all subprocesses will run
  sequentially.

## Value

A shiny app. This function returns a shiny app that can be run to
interact with the model.

## Details

If you set `nw = 3` and `nsw = 4`, you should have at least 16 cores,
one core for the shiny main process. Three cores for the three worker
processes and 12 cores (3 \* 4) for the subworkers. For the default
case, `nworkers = 2` and `nsw = 1`, you only need 3 cores, as `nsw = 1`
means that all subprocesses will run sequentially.

## Examples

``` r
if (interactive()) start_gui()
```
