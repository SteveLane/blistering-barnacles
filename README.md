# Risk factors for fouling biomass: Evidence from small vessels in Australia

This project contains the scripts to generate the paper.

## Reproducing analysis

All analyses and the accompanying manuscript can be easily reproduced using [docker](https://www.docker.com/) for a consistent build environment. The steps required to build the manuscript are:

1. Install and run docker.
2. Pull the required docker image (instructions forthcoming).
3. `make` the analysis magic happen!

### Installing and running docker

Head over to [https://www.docker.com/](https://www.docker.com/) to install docker for your operating system. Docker is a container platform, that enables users to run multiple different operating systems on top of their day-to-day os. This makes it easy to provide a consistent environment for data analysis: we can choose the specific os and package versions to install.

### Get the required docker image

Once you have docker up and running, pull the docker image into your system:

```sh
$ docker pull stevelane/analysis-dockerfiles:blistering-barnacles
```

(If you're not using the command line, you can use [kitematic](https://kitematic.com/) to search for the above image.)

This command pulls down the required software to reproduce the analysis. This consists of the operating system, R (and the various packages used) and latex.

### `make` the analysis

We're almost ready to `make` the analysis. First, we need to run the container that will do the work:

```sh
$ docker run -v $(pwd):/home/blistering-barnacles -it \
  stevelane/analysis-dockerfiles:blistering-barnacles
```

Here, `/path/on/your/machine` refers to where you would like to store the code and all the output from the analysis. When you run the above command, it will download the data from git repository and make it available for use.

Once the container has started up, you should be presented with a bash shell, in the appropriate starting directory. The image hash (the jumble of numbers and letters) will be different from that shown below, but your entrypoint will look something like:

```
root@6f5865b66643:/home/blistering-barnacles#
```

At this stage, you can just run `make`:

```
root@6f5865b66643:/home/blistering-barnacles# make
```

**NOTE**: by default, this will perform a multiple imputation analysis, using 50 multiply imputed datasets, and fit 5 MCMC models with 2000 iterations to each. It will take a long time! Run this on a machine which has lots of cores, or be prepared to wait. Alternatively, to test that it compiles, you can pass the number of multiply imputed datasets and number of iterations through as variables to the `make` script. E.g. the following will use 2 multiply imputed datasets, and 200 iterations:

```
root@6f5865b66643:/home/blistering-barnacles# make NUMMI=2 MCITER=200
```

--- **alternatively** ---

the graphics required to compile the manuscript are included as part of the repository. To compile the manuscript only (without running the models), use 

```
root@6f5865b66643:/home/blistering-barnacles# make paper
```

Let me know if you have any problems by filing an [issue](https://github.com/SteveLane/blistering-barnacles/issues).
