# SHE CTE

## Software identification
* Processing Element Name: PF-SHE
* Project Name: SHE_CTE (Common testing Environment)
* Profile: develop
* Version: 8.3 (10/11/2021)

## Contributors
### Active Contributors

* Bryan Gillis (b.gillis@roe.ac.uk)
* Rob Blake (rpb@roe.ac.uk)
* Gordon Gibb (gordon.gibb@ed.ac.uk)

### Other Contributors

* Giuseppe Congedo (giuseppe.congedo@ed.ac.uk)
* Niraj Welikala
* Nick Cross (njc@roe.ac.uk)
* Christoper Duncan (christopher.duncan@physics.ox.ac.uk)

## Purpose
This project contains various executables which constitute components in the shear estimation, validation and calibration pipelines.


## Relevant Documents
> `Boilerplate section which links any Euclid related documents that are relevant for the project`

* [RSD](link here)
* [SDD]()
* [VP/STS]()
* [STP/STR]()
* Any other relevant documents

## Dependencies

### Internal Euclid Dependencies

* [SHE_LensMC 3.3](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_LensMC)
* [SHE_MomentsML 8.2](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_MomentsML)
* [SHE_PPT 8.11 (indirect)](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT)


### External Euclid Dependencies

* [EL_Utils 1.1.0](https://gitlab.euclid-sgs.uk/EuclidLibs/EL_Utils)
* [ST_DataModelTools 8.0.5](https://gitlab.euclid-sgs.uk/ST-DM/ST_DataModelTools)
* [ST_DataModelBindings 8.0.5](https://gitlab.euclid-sgs.uk/ST-DM/ST_DataModelBindings)
* [ST_DataModel 8.0.5](https://gitlab.euclid-sgs.uk/ST-DM/ST_DataModel)
* [Elements 5.15](https://gitlab.euclid-sgs.uk/ST-TOOLS/Elements)


### Configuration

**EDEN 2.1**
```
- astropy 3.2.1
- numpy 1.17.2
- etc
```

### Dependant Projects

* [SHE_IAL_Pipelines](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines)


### Dependant Pipelines

* [SHE Analysis](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) 
* [Shear Calibration](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py)
* [SHE Global Validation](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Global_Validation/PipDef_SHE_Global_Validation.xml)
* [SHE_Scaling_Experiments](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py)
* [SHE_Shear_Reconciliation](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Reconciliation/PipScript_SHE_Shear_Reconciliation.py)



## Installation

All Euclid projects will be deployed via cvmfs. If this is installed and set up, this project will be pre-installed and no further work will be necessary. In case cvmfs isn't installed, or you wish to install an unreleased build or older build, you can do so through the following process:

```bash
cd ${HOME}/Work/Projects
git clone git@gitlab.euclid-sgs.uk:PF-SHE/SHE_CTE.git
cd SHE_CTE
git checkout <desired branch or tag>
make
make test
make install
```

## Main Programs Available

* SHE_CTE_BiasMeasurement: Executables for bias measurement
  - `SHE_CTE_MeasureStatistics`: For a set of shear measurements, computes statistics on them
  - `SHE_CTE_MeasureBias`: Given a set of statistics files (from `SHE_CTE_MeasureStatistics`), calculates the bias
  - `SHE_CTE_PlotBias`: Plots the bias as determined by `SHE_CTE_MeasureBias`
  - `SHE_CTE_PlotPsfSensitivity`: Plots the PSF sensitivity
  - `SHE_CTE_PrintBias`: Prints the bias as determined by `SHE_CTE_MeasureBias`
* SHE_CTE_PipelineUtility: Miscellaneous helper executables
  - `SHE_CTE_ObjectIdSplit`: Splits a list of objects from the input catalogue into a number of smaller batches
  - `SHE_CTE_SubObjectIdSplit`: Splits batches of objects into even smaller batches (see `SHE_CTE_ObjectIdSplit`)
  - `SHE_CTE_ShearEstimatesMerge`: ??
  - `SHE_CTE_CleanupBiasMeasurement`: Cleans up intermediate files involved in bias measurement
* SHE_CTE_ScalingExperiments: Executables for a test pipeline to determine how to improve the scaling of the analysis pipeline
  - `SHE_CTE_SplitFits`: splits fits file for an exposure into individual FITS/HDF5 files, one per CCD
  - `SHE_CTE_CombineSplitFitsListfile`: Combines output from `SHE_CTE_SplitFits` into a single json
  - `SHE_CTE_ExtractObjects`: For an input catalogue and observation, extracts all the objects from the catalogue that are present in the observation
  - `SHE_CTE_MakeBatches`: Splits the list of objects in the observation into small batches
  - `SHE_CTE_ExtractStamps`: Given a list of objects, times itself extracting 400x400 pixel postage stamps around the objects
  - `SHE_CTE_AnalyseRuntime`: Analyses the runtime of the pipeline, producing graphs
* SHE_CTE_ShearEstimation: Shear Estimation executables
  - `SHE_CTE_EstimateShear`: Given a list of objects, estimates their shear
* SHE_CTE_ShearReconcilliation: Shear Reconciliation Executables
  - `SHE_CTE_ReconcileShear`: ??
* SHE_CTE_ShearValidation: Shear validation executables
  - `SHE_CTE_CrossValidateShear`: ??


## Running the software

### `SHE_CTE_MeasureStatistics`

This executable takes the actual and measured shears for each object and from these constructs and stores shear statistics data tables for them (see the LinregressStatistics object in [`SHE_PPT/math.py`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/math.py))

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_MeasureStatistics --details_table <file> --shear_estimates <file> --pipeline_config <file> --she_bias_statistics <file> [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --details_table `<filename>`     | `.xml` data product pointing to the details table (OU-MER Final Catalogue??) | yes    | N/A         |
| --shear_estimates `<filename>`  | `.xml` data product pointing to the shear measurement tables for each measurement method   | yes | N/A |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --she_bias_statistics `<filename>`      | `.xml` data product pointing to bias statistics tables. | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`details_table`_:

**Description:** The filename of an `.xml` MerFinalCatalog data product in the workdir. This data product and the format of the associated data table are described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html.

**Source:**  Generated by the `SHE_GST_GenGalaxyImages` executable from [`SHE_GST`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST)


_`shear_estimates`_:

**Description** The filename of a sheMeasurements `.xml` data product in the workdir. This data product and the format of the associated data table are described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html.

**Source** Generated by `SHE_CTE_EstimateShear`

_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default behaviour** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_MeasureStatistics_archive_dir `<dir>` | The directory for archiving the data.   | N/A |
| SHE_CTE_MeasureStatistics_webdav_archive `<True/False>`  | Sets whether to archive to webdav or not | Does not archive to webdav |
| SHE_CTE_MeasureStatistics_webdav_dir `<dir>` | The webdav directory to archive to | N/A |

If both these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_bias_statistics`_ 

**Description** The desired name of the output shear statistics `.xml` data product

**Details** The product is of the type sheIntermediateGeneral, and points to a shear bias statistics table FITS file (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py))


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.













### `SHE_CTE_MeasureBias`

From a set of input bias statistics products, combines these and calculates the bias for each measurement method. Also has a "recovery mode" whereby it will search the workdir for statistics products to use.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_MeasureBias --she_bias_statistics <file>  --pipeline_config <file> --she_bias_measurements <file> [--details_table_head <str>] [--bootstrap-seed <int>] [--recovery_bias_statistics_filename <filename>] [--recovery_bias_measurements_filename <filename>] [--store_measurements_only] [--use_bias_only] [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --she_bias_statistics `<filename>`     | Listfile pointing to a number of bias statistics data products | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --she_bias_measurements `<filename>`  | Name of the `.xml` sheBiasStatistics data product containing the bias measurements | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --bootstrap_seed `<int: seed>`    | The random seed for bootstrapping the bias measurement errors  | no           | 0       |
| --details_table_head `<str>`    | "Desired head for the filenames of the output details tables" (NOT USED IN CODE AT ALL)  | no           | None       |
| --recovery-bias-statistics-filename `<filename>`    | The filename to look for when looking for bias statistics files in recovery mode  | no           | `she_bias_statistics.xml`       |
| --recovery-bias-measurements-filename `<filename>`    | The filename to look for when looking for bias measurements files in recovery mode  | no           | `she_bias_measurements.xml`       |
| --store_measurements_only ("store_true")    | Only store the bias measurements, and none of the intermediary bias statistics data   | no           | False      |
| --use_bias_only ("store_true")    | Only use existing bias measurements to calculate the bias, ignoring any bias statistics   | no           | False      |
| --number_threads `<int: nthreads>`    | Number of "threads" (actually processes) to use when reading in the input files   | no           | 8      |
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`she_bias_statistics`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` sheIntermediateGeneral data products, that point to a shear bias statistics FITS tables (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py))

**Source:**  Files pointed to in the listfile are generated by the `SHE_CTE_MeasureStatistics` executable. The listfile is generated automatically by the shear calibration pipeline when combining parallel branches of the pipeline together.


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default behaviour** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_MeasureBias_archive_dir `<dir>` | The directory for archiving the data.   | N/A |
| SHE_CTE_MeasureBias_webdav_archive `<True/False>`  | Sets whether to archive to webdav or not | Does not archive to webdav |
| SHE_CTE_MeasureBias_webdav_dir `<dir>` | The webdav directory to archive to | N/A |

If both these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_bias_measurements`_ 

**Description** The desired name of the output bias measurements `.xml` data product

**Details** The product is of the type sheIntermediateGeneral, and points to a shear bias statistics table FITS file (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py)) which contains the bias measurements in its header, but also can contain the bias statistics in the body of the table. 


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_PlotBias`

Plots Bias? To do with sensitivity? This executable does not seem to be run in any pipeline, but is instead run manually in post-processing.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_PlotBias --methods <method list>  --root-data-folder <dir> --bias_measurements_head <str> --output_file_name_head <str> [--output_format <str>] [--hide] [--plot_error] [--plot_slopes] [--normed_only] [--unnormed_only] [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --methods `<str>`     | A list of methods that we want to generate plots for | no    |  ["KSB", "REGAUSS", "MomentsML", "LensMC"]       |
| --root_data_folder `<dir>` | A directory that contains bias measurements | no          | same as workdir         |
| --bias_measurements_head `<str>` | Head of filenames of bias measurements, minus the E??P??S?? portion | no          | `shear_bias_measurements`         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_file_name_head `<str>` | Head of output filenames  | no          | `sensitivity_testing`         |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_format `<str>` | Output plot file format | no          | `png`         |
| --hide (`store true`)    | If set, will not display plots to the screen  | no           | false       |
| --plot_error (`store true`)    | If set, will plot the errors for all parameters  | no           | false       |
| --plot_slopes (`store true`)    | If set, will plot the slopes and slope errors for all methods  | no           | false       |
| --normed_only (`store true`)    | If set, only show normed plots  | no           | false       |
| --unnormed_only (`store true`)    | If set, will only show unnormed plots  | no           | false       |
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |



_**Outputs**_

Presumably a series of plots of various quantities, but this is not documented anywhere...


_**Example**_

Unknown








### `SHE_CTE_PlotPsfSensitivity`

Plots PSF Sensitivity? This executable does not seem to be run in any pipeline, but is instead run manually in post-processing.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_PlotPsfSensitivity` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_PlotPsfSensitivity --methods <method list>  --root-data-folder <dir> --bias_measurements_head <str> --output_file_name_head <str> [--output_format <str>] [--hide] [--plot_error] [--normed_only] [--unnormed_only] [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --methods `<str>`     | A list of methods that we want to generate plots for | no    |  ["KSB", "REGAUSS", "MomentsML", "LensMC"]       |
| --root_data_folder `<dir>` | A directory that contains bias measurements | no          | same as workdir         |
| --bias_measurements_head `<str>` | Head of filenames of bias measurements, minus the E??P??S?? portion | no          | `shear_bias_measurements`         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_file_name_head `<str>` | Head of output filenames  | no          | `sensitivity_testing`         |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_format `<str>` | Output plot file format | no          | `png`         |
| --hide (`store true`)    | If set, will not display plots to the screen  | no           | false       |
| --plot_error (`store true`)    | If set, will plot the errors for all parameters  | no           | false       |
| --normed_only (`store true`)    | If set, only show normed plots  | no           | false       |
| --unnormed_only (`store true`)    | If set, will only show unnormed plots  | no           | false       |
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |



_**Outputs**_

Presumably a series of plots of various quantities, but this is not documented anywhere...


_**Example**_

Unknown



### `SHE_CTE_PrintBias`

Prints the bias calculated by `SHE_CTE_MeasureBias` out to screen.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_PrintBias` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_PlotPsfSensitivity --she_bias_measurements <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --she_bias_measurements `<file>` | Filename of the bias measurements file | yes          | N/A         |


**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |



_**Inputs**_

_`she_bias_measurements`_:

**Description:** The filename of an `.xml` sheIntermediateGeneral data product, that points to shear bias statistics FITS tables (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py)) that contain bias measurements in their headers

**Source:**  The bias measurement file is produced by `SHE_CTE_MeasureBias` executable during execution of the Calibration pipeline.


_**Example**_

Unknown









### `SHE_CTE_CleanupBiasMeasurement`

Cleans up intermediate files that are created during the pipeline execution, but are not needed any more.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_CleanupBiasMeasurement` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_CleanupBiasMeasurement --simulation_config <file> --data_images <file> --stacked_data_image <file> --psf_images_and_tables <file> --segmentation_images <file> --stacked_segmentation_image <file> --detections_tables <file> --details_table <file> --shear_estimates <file> --shear_bias_statistics_in <file> --shear_bias_statistics_out <file> --pipeline_config <file>  [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --simulation_config `<filename>`     | XML data product containing the simulation configuration | yes    | N/A         |
| --data_images `<filename>`     | Listfile pointing to XML data products for the individual exposures | yes    | N/A         |
| --stacked_data_image `<filename>`     | XML data product for the stacked image | yes    | N/A         |
| --psf_images_and_tables `<filename>`     | Listfile pointing to PSF data products | yes    | N/A         |
| --segmentation_images `<filename>`     | Listfile pointing to XML segmentation map data products for the exposures | yes    | N/A         |
| --stacked_segmentation_image `<filename>` | XML data product for the stacked image's segmentation map | yes    | N/A         |
| --detections_tables `<filename>` | Listfile pointing to a number of detection table XML data products | yes    | N/A         |
| --details_table `<filename>` | XML data product for the details table  | yes    | N/A         |
| --shear_estimates `<filename>` | XML data product for the shear measurements | yes    | N/A         |
| --shear_bias_statistics_in `<filename>` | XML data product for the bias statistics | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --shear_bias_statistics_out `<filename>` | Name for the output bias statistics product (NB this is just a copy of `shear_bias_statistics_in`) | yes    | N/A         |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`simulation_config`_:

**Description:** The filename of a `.xml` data product containing configuration options for `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))

**Source:**  Created by `SHE_GST_PrepareConfigs` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`data_images`_:

**Description:** A listfile containing [visCalibratedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_calibratedframe.html) data products.

**Source:**  These are generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`stacked_data_image`_:

**Description:** An `.xml` [visStackedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_visstackedframe.html) data product.

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`psf_images_and_tables`_:

**Description:** An `.xml` [shePsfModelImage](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_psfmodelimage.html) data product.

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`segmentation_images`_:

**Description:** A listfile pointing to `.xml` [sheExposureReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_exposurereprojectedsegmentationmap.html) data products for each exposure.

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`stacked_segmentation_image`_:

**Description:** An `.xml` [sheStackReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_stackreprojectedsegmentationmap.html) data product.

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`detections_tables`_:

**Description:** A Listfile pointing to detection table data products?

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`details_table`_:

**Description:** An `.xml` details_table data product?

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`shear_estimates`_:

**Description:** An `.xml` [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html) data product.

**Source:**  Generated by `SHE_CTE_EstimateShear`


_`shear_bias_statistics_in`_:

**Description:** An `.xml` sheIntermediateGeneral, and points to a shear bias statistics table FITS file (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py))

**Source:**  Generated by `SHE_CTE_MeasureStatistics`


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default behaviour** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_CleanupBiasMeasurement_cleanup | Option to stop this executable removing the files (if set to False)   | True |


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_bias_statistics_out`_ 

**Description** The desired name of the output bias measurements `.xml` data product. This is a copy of `she_bias_statistics_in`, required so that this executable outputs something and hence the pipeline runner will actually run this executable in the pipeline.

**Details** The product is of the type sheIntermediateGeneral, and points to a shear bias statistics table FITS file (described [here](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py)).


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data, *_although as its job is to delete files it may throw an exception as the files it has to delete will no longer exist!_* The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_ObjectIdSplit`

Finds objects in the MER catalogue that are in the observation, then divides them into batches for processing by other executables. This is part of the analysis pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ObjectIdSplit` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_ObjectIdSplit --mer_final_catalog_tables <file> --data_images <file>  --pipeline_config <file> --object_ids <file> --batch_mer_catalogs <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --mer_final_catalog_tables `<filename>`     | Listfile pointing to a number of MER final catalogues | yes    | N/A         |
| --data_images `<filename>`     | Listfile pointing to a number of VIS exposures | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --object_ids `<filename>`  | Name of the `.json` listfile containing the list of sheObjectIdList data products for each batch | yes          | N/A |
| --batch_mer_catalogs `<filename>`  | Name of the `.json` listfile containing the list of merFinalCatalog products for each batch | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`mer_final_catalog_tables`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) data products.

**Source:** This is provided as input to the Analysis pipeline


_`data_images`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [visCalibratedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_calibratedframe.html) data products.

**Source:** This is provided as input to the Analysis pipeline


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_ObjectIdSplit_batch_size `<int>` | The number of objects per batch.   | 400 |
| SHE_CTE_ObjectIdSplit_max_batches `<int>`  | The maximum number of batches to use. If it's 0, the number of batches is unlimited | 0 |
| SHE_CTE_ObjectIdSplit_batch_ids `<str?>` | A list of object IDs to use | None |


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`object_ids`_ 

**Description** The desired name for a listfile pointing to numerous [sheObjectIdList](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_objectidlist.html) products, one per batch.

_`batch_mer_catalogs`_ 

**Description** The desired name for a listfile pointing to numerous [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) products, one per batch. These contain the same information as the input MER final catalogues, just only the rows for each object in the batch.

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.





### `SHE_CTE_SubObjectIdSplit`

This executable is virtually identical to SHE_CTE_ObjectIdSplit, save for its pipeline_config options.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_SubObjectIdSplit` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_SubObjectIdSplit --mer_final_catalog_tables <file> --data_images <file>  --pipeline_config <file> --object_ids <file> --batch_mer_catalogs <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --mer_final_catalog_tables `<filename>`     | Listfile pointing to a number of MER final catalogues | yes    | N/A         |
| --data_images `<filename>`     | Listfile pointing to a number of VIS exposures | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --object_ids `<filename>`  | Name of the `.json` listfile containing the list of sheObjectIdList data products for each batch | yes          | N/A |
| --batch_mer_catalogs `<filename>`  | Name of the `.json` listfile containing the list of merFinalCatalog products for each batch | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`mer_final_catalog_tables`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) data products.

**Source:** This is provided as input to the Analysis pipeline


_`data_images`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [visCalibratedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_calibratedframe.html) data products.

**Source:** This is provided as input to the Analysis pipeline


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_SubObjectIdSplit_batch_size `<int>` | The number of objects per batch.   | 20 |
| SHE_CTE_SubObjectIdSplit_max_batches `<int>`  | The maximum number of batches to use. If it's 0, the number of batches is unlimited | 0 |
| SHE_CTE_SubObjectIdSplit_batch_ids `<str?>` | A list of object IDs to use | None |


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`object_ids`_ 

**Description** The desired name for a listfile pointing to numerous [sheObjectIdList](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_objectidlist.html) products, one per batch.

_`batch_mer_catalogs`_ 

**Description** The desired name for a listfile pointing to numerous [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) products, one per batch. These contain the same information as the input MER final catalogues, just only the rows for each object in the batch.

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_ShearEstimatesMerge`

Given input lists shear estimates and LensMC chains, combines them into two tables, one for the estimates, and one for the chains. This is a merge point executable in the Analysis pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ShearEstimatesMerge` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_ShearEstimatesMerge --shear_estimates_product_listfile <file> --she_lensmc_chains_listfile <file>  --pipeline_config <file> ----merged_she_measurements <file> --merged_she_lensmc_chains <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --shear_estimates_product_listfile `<filename>`     | Listfile pointing to a number of shear measurements | yes    | N/A         |
| --she_lensmc_chains_listfile `<filename>`     | Listfile pointing to a number of LensMC chains products | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --merged_she_measurements `<filename>`  | Name of the `.xml` shear measurement product | yes          | N/A |
| --merged_she_lensmc_chains `<filename>`  | Name of the `.xml` LensMC chains product | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`shear_estimates_product_listfile`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html) data products.

**Source:** This is produced by SHE_CTE_EstimateShear in the Analysis pipeline


_`she_lensmc_chains_listfile`_:

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` [sheLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmcchains.html) data products.

**Source:** This is produced by SHE_CTE_EstimateShear in the Analysis pipeline


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_SHE_CTE_ShearEstimatesMerge_number_threads `<int>` | The number of "threads" (processes) to use when reading in the input files.   | 8 |


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`merged_she_measurements`_ 

**Description** The desired name for the merged [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html) data product.

_`merged_she_lensmc_chains`_ 

**Description** The desired name for the merged [sheLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmcchains.html) data product.

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_SplitFits`

Takes the input exposure and splits it into individual files, one per CCD. This executable can produce HDF5 or FITS files depending upon the configuration options set. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_SplitFits` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_SplitFits --input_fits_json <file> --pipeline_config <file> --output_json <file> --timing_info <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --input_fits_json `<filename>`     | json that points to the individual fits files for that exposure | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_json `<filename>`  | Name of the output json, which points to the new CCD files | yes          | N/A |
| --timing_info `<filename>`  | Name of the json file to write timing statistics to | yes          | N/A |



_**Inputs**_

_`input_fits_json`_:

**Description:** The filename of a `.json` file pointing to the VIS FITS files for an exposure. An example of the file would be:
```
{
    "exp": "1",
    "det": "/path/to/detector/fits/file/for/exposure.fits",
    "wgt": "/path/to/weight/fits/file/for/exposure.fits",
    "bkg": "/path/to/background/fits/file/for/exposure.fits"
}
```

**Source:** This is produced manually as input.



_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| HDF5 `<int>` | Create HDF5 files if this flag is not equal to zero, otherwise create FITS files   | 1 |
| chunked `<int>` | Create chunked HDF5 files if this flag is not equal to zero, otherwise do not use chunking. This option only comes into play if we are creating HDF5 files   | 1 |
| memmap `<int>` | Open input FITS files with memmap   | 1 |
| compression `<int>` | Use compression when writing HDF5 files This option only comes into play if we are creating HDF5 files   | 0 |

**Source:** 

User Generated

_**Outputs**_

_`output_json`_ 

**Description** The desired name for the output json file that points to the newly created FITS/HDF5 files for the CCDs. The file is of the form:
```
{
    "1": {
        "1-1": "/path/to/file/for/ccd_1-1",
        "1-2": "/path/to/file/for/ccd_1-2",
        "1-3": "/path/to/file/for/ccd_1-3",
        ...,
        "6-6": "/path/to/file/for/ccd_6-6",
    }
}
```

_`timing_info`_ 

**Description** The desired name for the output timing information json. This file is of the form:
```
{
    "tstart": "YYYY-MM-DD HH:MM:SS.xxxxxx",
    "walltime": 100.0,
    "exposure": 1
}
```

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.




### `SHE_CTE_CombineSplitFitsListfile`

Combines the output from SHE_CTE_SplitFits into a single json. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_CombineSplitFitsListfile` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_CombineSplitFitsListfile --input_listfile <file> --pipeline_config <file> --output_json <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --input_listfile `<filename>`     | Listfile pointing to the output_json files from the runs of SHE_CTE_SplitFits  | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_json `<filename>`  | Name of the output json combining the inputs | yes          | N/A |




_**Inputs**_

_`input_listfile`_:

**Description:** The filename of a `.json` listfile pointing to the outputs from SHE_CTE_SplitFits runs, one per exposure


**Source:** This is produced automatically by the pipeline runner from individual executions of SHE_CTE_SplitFits



_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs. At present there are no options in this file that are used by SHE_CTE_CombineSplitFitsListfile

**Source:** 

User Generated



_**Outputs**_

_`output_json`_ 

**Description** The desired name for the output json file contains the combined output from the SHE_CTE_SplitFits steps. The file is of the form:
```
{
    "1": {
        "1-1": "/path/to/file/for/exp1/ccd_1-1",
        "1-2": "/path/to/file/for/exp1/ccd_1-2",
        "1-3": "/path/to/file/for/exp1/ccd_1-3",
        ...,
        "6-6": "/path/to/file/for/exp1/ccd_6-6",
    },
    "2": {
        "1-1": "/path/to/file/for/exp2/ccd_1-1",
        "1-2": "/path/to/file/for/exp2/ccd_1-2",
        "1-3": "/path/to/file/for/exp2/ccd_1-3",
        ...,
        "6-6": "/path/to/file/for/exp2/ccd_6-6",
    },
    ...
}
```


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.



### `SHE_CTE_ExtractObjects`

Given an input catalogue of objects and an observation, extracts the objects from the catalogue that are in the observation, recording which CCD the object is on for each exposure. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ExtractObjects` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_ExtractObjects --stacked_image <file> --exposures <file> --catalogue_listfile <file> --pipeline_config <file> --output_objects_list <file> --combined_catalogue <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --stacked_image `<filename>`     | FITS file for the VIS stacked image for the observation  | yes    | N/A         |
| --exposures `<filename>`     | JSON file pointing to the files for each CCD for each exposure  | yes    | N/A         |
| --catalogue_listfile `<filename>`     | Listfile for the MER catalogue FITS files  | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_objects_list `<filename>`  | JSON containing a list of all the objects in the observation | yes          | N/A |
| --combined_catalogue `<filename>`  | FITS table of all the input MER catalogues added together | yes          | N/A |



_**Inputs**_

_`stacked_image`_:

**Description:** The filename of the VIS stacked image FITS file for an observation

**Source:** This is supplied as an input to the pipeline

_`exposures`_:

**Description:** JSON pointing to the individual files for each CCD for each exposure. 

**Source:** This is the output from SHE_CTE_CombineSplitFitsListfile

_`catalogue_listfile`_:

**Description:** JSON pointing to a number of MER Catalogue FITS files. 

**Source:** This is supplied as input to the pipeline

_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs. At present there are no options in this file that are used by SHE_CTE_ExtractObjects

**Source:** 

User Generated



_**Outputs**_

_`output_objects_list`_ 

**Description** The desired name for the output json file containing a list of the objects found in the observation, along with some of their properties. The file is of the form:
```
[
  {
        "1": "4-3",
        "2": "4-3",
        "3": "4-3",
        "4": "4-3",
        "ra": 232.14187741659418,
        "dec": 30.237904069815297,
        "index": 215811,
        "x": 16525,
        "y": 22285
  },
  {
        "1": "4-1",
        "3": "5-1",
        "4": "5-1",
        "ra": 231.88449735119158,
        "dec": 30.240277551936423,
        "index": 215857,
        "x": 9044,
        "y": 25061
  },
  ...
]
```
With the properties being the detector the object is in for each exposure, the object's RA, Dec, its index in the output table (see below) and the x and y pixels of the object on the stacked image.

_`combined_catalogue`_ 

**Description** The desired name for the FITS table combining all the input catalogues into one.

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.




### `SHE_CTE_MakeBatches`

Takes the list of all the objects in the observation and splits them into small batches. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MakeBatches` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_MakeBatches --objects_list <file> --pipeline_config <file> --batch_listfile <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --objects_list `<filename>`     | json that contains a list of all the objects in the exposure | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --batch_listfile `<filename>`  | Name of the output listfile, which points to a number of new objects list files, one per batch | yes          | N/A |



_**Inputs**_

_`objects_list`_:

**Description:** The filename of a `.json` file containing the list of all the objects in the observation.


**Source:** This is produced by SHE_CET_ExtractObjects



_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| maxbatches `<int>` | Maximum number of batches to create (-1 means unlimited number of batches)   | 4 |
| batchsize `<int>` | The (mean) number of objects per batch   | 20 |
| memmap `<int>` | Open input FITS files with memmap   | 1 |
| spatial_batching `<int>` | Whether to spatially batch the objects (1) or batch them by the order they come in the input objects list (0)   | 0 |

**Source:** 

User Generated

_**Outputs**_

_`batch_listfile`_ 

**Description** The desired name for the output json file that points to the newly created object lists for each batch. The object lists are the same format as the `output_objects_list` from SHE_CTE_ExtractObjects.


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.







### `SHE_CTE_ExtractStamps`

Given the input list of objects, extracts 400x400 pixel postage stamps from each CCD the object is in and times how long it takes to do this. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ExtractStamps` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_ExtractStamps --batch_file <file> --exposures <file> --pipeline_config <file> --timing_info <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --batch_file `<filename>`     | json that contains a list of all the objects to extract stamps for | yes    | N/A         |
| --exposures `<filename>`     | json that points to all the FITS/HDF5 files for the individual CCDs | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --timing_info `<filename>`  | json file that contains timing statistics for the execution of this executable | yes          | N/A |



_**Inputs**_

_`batch_file`_:

**Description:** The filename of a `.json` file containing the list of all the objects to be read in. This is the same format as `output_objects_list` from SHE_CTE_ExtractObjects


**Source:** This is produced by SHE_CTE_MakeBatches


_`exposures`_:

**Description:** The filename of the `.json` file containing the list of all the FITS/HDF5 files for the individual CCDs for each exposure. This is the same format as `output_json` from SHE_CTE_CombineSplitFitsListfile


**Source:** This is produced by SHE_CTE_CombineSplitFitsListfile


_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default value** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| mean_compute_time `<int>` | Mean time in seconds to pause for after reading a stamp in to simulate the time taken to measure the shear for that object   | 10 |


**Source:** 

User Generated

_**Outputs**_

_`timing_info`_ 

**Description** The desired name for the output json file that contains timing information for the execution. The file is of the form:
```
{
    "num_objects": 21,
    "num_files": 4,
    "walltime": 13.528909921646118,
    "tstart": "2021-09-09 10:52:27.536006",
    "compute_time": 0.0
}
```
Where it provides the number of objects that had stamps extracted, the number of files (CCDs) that were used in reading in the objects, the walltime of the executable, the start time of the executable and the time taken simulating compute.


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.







### `SHE_CTE_AnalyseRuntime`

Analyses the runtime of SHE_CTE_SplitFITS and SHE_CTE_ExtractStamps and produces a number of graphs and a json file with some statistics. This is the final step of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_AnalyseRuntime` on Elements use the following command:
```bash
E-Run SHE_CTE 8.3 SHE_CTE_AnalyseRuntime --stamp_timing_listfile <file> --split_timing_listfile <file> --pipeline_config <file> --results <file> [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |
| --log-dir `<file>`    | The name of the file to output the executable logs to. If not provided, logging data will only be output to the terminal. When run via the pipeline runner, this will be set to a file in the directory `<workdir>/logs/<task_name>/` with a name based off of the command used to call this executable, such as "E_Run_SHE_MyProject_0.1_SHE_MyProject_GenCatPic_retry_0.out" | no           | .        |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --stamp_timing_listfile `<filename>`     | json that contains a lisfile of all the timing_info outputs from SHE_CTE_ExtractStamps | yes    | N/A         |
| --split_timing_listfile `<filename>`     | json that contains a lisfile of all the timing_info outputs from SHE_CTE_SplitFits | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --results `<filename>`  | json file that contains timing statistics for the execution of this executable and points to plots produced by the executable | yes          | N/A |



_**Inputs**_

_`stamp_timing_listfile`_:

**Description:** The filename of a `.json` file containing the list of all the `timing_info` outputs from SHE_CTE_ExtractStamps.


**Source:** This is produced by the pipeline runner from the outputs from SHE_CTE_ExtractStamps

_`split_timing_listfile`_:

**Description:** The filename of a `.json` file containing the list of all the `timing_info` outputs from SHE_CTE_SplitFits.


**Source:** This is produced by the pipeline runner from the outputs from SHE_CTE_SplitFits


_`pipeline_config`_:

**Description:**

A plain text file containing key-value pairs.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables. None of the options in this file are currently used by the executable

**Source:** 

User Generated


_**Outputs**_

_`results`_ 

**Description** The desired name for the output json file that contains timing statistics for the pipeline execution, and links to the plots profuced. The file is of the form:
```
{
    "Total runtime for whole pipeline": 425.227296,
    "Total SplitFits Walltime": 49.106916,
    "Mean SplitFits task walltime": 0.02664393186569214,
    "Total SplitFits CPUh": 2.960436873965793e-05,
    "Total ExtractStamps walltime": 137.093798,
    "Mean walltime for ExtractStamps tasks": 15.826434627175331,
    "Mean ExtractStamps task I/O time": 15.826434627175331,
    "Min ExtractStamps task I/O time": 7.03668999671936,
    "Max ExtractStamps task I/O time": 23.88964295387268,
    "Standard Deviation of ExtractStamps task I/O time": 4.745512086551991,
    "ExtractStamps I/O CPUh": 0.14067941890822516,
    "Total ExtractStamps CPUh": 0.14067941890822516,
    "Mean number of running ExtractStamps tasks at any given time": 3.7007299270072993,
    "Max number of running ExtractStamps tasks at any given time": 7,
    "I/O time Vs Num Objects": "/path/to/iotime_objects.png",
    "I/O time Vs Num Files": "/path/to/iotime_files.png",
    "Start Time Vs I/O time": "/path/to/start_time_stamp.png",
    "Total duration Vs Batch Number": "/path/to/duration_batchno.png",
    "Number of running tasks with time": "/path/to/num_running_tasks.png",
    "Total duration Vs Exposure": "/path/to/duration_split.png",
    "split_stamp_walltimess for Exposures": "/path/to/split_stamp_exposure_walltimes.png",
    "pipeline_config": {
        "HDF5": 0,
        "chunked": 1,
        "compression": 0,
        "maxbatches": 32,
        "batchsize": 20,
        "memmap": 1,
        "spatial_batching": 1,
        "mean_compute_time": 0,
        "dry_run": 0
    }
}
```

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.










### `SHE_MyProject_ShowCatPic`
> `Same structure as before: how to run the code on Elements, what are the options for the command line with descriptions and what each external file and a simple example for the user to run`

## Troubleshooting
> `If any problems are anticipated, add a section here for them to help users resolve them on their own before they need to appeal to a developer for help.`

### The cat in the generated picture appears to be wearing both a standard tie and a bowtie
> `For known problems which can occur if the user makes a common error, explain how it can be resolved.`

This is a known bug which occurs if the user requests `--use_tie=bowtie`. The correct argument is `--use_tie=bow`.

### A test failed when I ran "make test"

_**Ensure you have the most up-to-date version of the project and all its dependencies**_

It's possible the issue you're hitting is a bug that's already been fixed, or could be due to locally-installed versions of projects on the develop branch no longer being compatible with a newly-deployed version of another dependency on CODEEN. If you're running on the develop branch and have installed locally, pull the project, call `make purge`, and install again, and repeat for all dependencies you've installed locally. Try running `make test` again to see if it works.

_**Report the failing test to the developers**_

If the test still fails, please report it to the active developers listed above, ideally by creating a GitLab issue, or else by e-mailing them.

_**Try running the desired code**_

Tests can fail for many reasons, and a common reason is that the code is updated but not the test. This could lead to the test failing even if the code works properly. After you've reported the issue, you can try to run the desired code before the issue with the failing test has been fixed. There's a decent chance that the bug might only be in the test code, and the executable code will still function.

### An exception was raised, what do I do?
> `General instructions for how to figure out what to do when an exception is raised, which can be copied for all projects.`

_**Check for an issue with the input**_

First, look through the exception text to see if it indicates an issue with the input data. This will often be indicated by the final raised exception indicating an issue with reading a file, such as a SheFileReadError which states it cannot open a file. If this is the case, check if the file exists and is in the format that the code expects. If the file doesn't exist, then you've found the problem. Either a needed input file is missing, or one of the input files points to the incorrect filename. Determine which this is, and fix it from there.

If the file does exist but you still see an error reading from it, then the issue is most likely that the file is unreadable for some reason - perhaps the download was corrupt, perhaps manual editing left it improperly formatted, etc. Try to test if this is the case by reading it manually. For instance, if the program can't open a `FITS` file, try opening it with `astropy`, `ds9`, `topcat` etc. (whatever you're comfortable with) to see if you can read it external to the code.

Keep in mind that the code might try multiple methods to open a file. For instance, the pipeline_config input file can be supplied as either a raw text file, an `.xml` data product, or a `.json` listfile. The program will try all these options, and if all fail, the final exception text will only show the final type attempted. The full traceback, however, should show all attempts. So if it appears that the program tried to read a file as the wrong type, check through the traceback to see if it previously tried to read it as the expected type and failed.

_**Ensure you have the most up-to-date version of the project and all its dependencies**_

It's possible the issue you're hitting is a bug that's already been fixed, or could be due to locally-installed versions of projects on the develop branch no longer being compatible with a newly-deployed version of another dependency on CODEEN. If you're running on the develop branch and have installed locally, pull the project, call `make purge`, and install again, and repeat for all dependencies you've installed locally. Try running again to see if this works.

_**See if the exception, traceback, or log gives you any other clue to solve the problem**_

There are many reasons something might go wrong, and many have been anticipated in the code with an exception to indicate this. The exception text might tell you explicitly what the problem is - for instance, maybe two options you set aren't compatible together. If it wasn't an anticipated problem, the exception text probably won't obviously indicate the source of the problem, but you might be able to intuit it from the traceback. Look through the traceback at least a few steps back to see if anything jumps out at you as a potential problem that you can fix. Also check the logging of the program for any errors or warnings, and consider if those might be related to your problem.

_**Report the issue**_

If all else fails, raise an issue with the developers on GitLab. Be sure to include the following information:

1. Any details of input data you're using.
1. The command you called to trigger the program (or the pipeline which called the program)
1. The full log of the execution, from the start of the program to the ultimate failure. In the case of a failure during a pipeline run, you can attach the generated log file for this executable, which can be found in the `logs` directory within the work directory, and then in a subdirectory corresponding to this task.
1. Any steps you've taken to try to resolve this problem on your own.

