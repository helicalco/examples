## Mrode examples

There are 4 examples from Mrode's textbook recreated in full using bash and helical tools demonstrating different basic models.

1. `03_1_univariate` - a univariate model
2. `04_1_repeatability` - a repeatability model
3. `05_1_multivariate` - a multivariate model for weaning weight and post weaning gain
4. `05_5_multivariate` - a multivariate model for fat 

Each of these scripts can be run independently to generate the complete set up and solution for these examples from Mrode.

Once any script has been run, you can `cd` into the associated run directory, e.g. `cd 03_1.run` for the `03_1_univariate` example to examine the outputs.

For each model we also include an example `helical bolt model` parameter file, e.g. `03_1.helical` for the `03_1_univariate` example. The parameter file requires the associated example to have been run, but can then be used to generate a complete set of turnkey model scripts to run the same exact model from scratch using only the raw data inputs and the parameter file.

To generate the turn-key scripts, first `cd 03_1.run` then run `helical bolt model ../03_1.helical`. A `scripts` folder will be created inside the run directory with a `runner` script that will run the analysis pipeline.  The `helical bolt model` tool outputs include an additional Gibbs sampling step using the BOLT `ssgibbsCuda` tool. In the near future Gibbs sampling will be available in the helical toolset and the `helical bolt model` tool will be updated to utilise this instead, however, it will leave drop-in replacement examples of how to use the BOLT `pcg` and `ssgibbsCuda` commands in place of the `helical` equivalent tools. The BOLT tools take advantage of GPU computing, while the `helical` toolset is CPU only.

Please contact Dan Garrick (dan@helicalco.com) with any questions or feedback, or please open an issue at our GitHub repository https://github.com/helicalco/examples
