# Graph Mining Notebooks

Julia Notebooks and datasets to accompany the textbook "Mining Networks" (working title).

## Software environment

Please make sure to install the Julia programming language and make its
exacutable discoverable from the command line.

Then, after cloning this repository to your local machine go to command line
and in the folder containing the Julia notebooks and Project.toml and Manifest.toml file
run the following command:

```
julia --project=. --eval "using Pkg; Pkg.instantiate()"
```

This will download and install all the packages required to work with our examples.

You can start Jupyter Notebooks either using `jupyter notebook` command
in the command line, or by running the:
```
julia --project=. --eval "using IJulia; notebook(dir=pwd())"
```
command (that starts Jupyter from within a Julia session).

Note that some of the notebooks require installing additional external libraries.
In this case please follow the additional installation instructions provided in them.
