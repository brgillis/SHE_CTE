# SHE CTE

## Software identification
* Processing Element Name: PF-SHE
* Project Name: SHE_CTE (Common Testing Environment)
* Profile: develop
* Version: 9.3 (2024/03/04)

## Contributors
### Active Contributors

* Bryan Gillis (b.gillis@roe.ac.uk)
* Rob Blake (rpb@roe.ac.uk)
* Gordon Gibb (gordon.gibb@ed.ac.uk)
* Richard Rollins (@rrollins)

### Other Contributors

* Giuseppe Congedo (giuseppe.congedo@ed.ac.uk)
* Niraj Welikala
* Nick Cross (njc@roe.ac.uk)
* Christoper Duncan (christopher.duncan@physics.ox.ac.uk)

## Purpose
This project contains various executables which constitute components in various SHE pipelines, in particular the analysis, calibration and validation pipelines. For info on these pipelines, please see the [SHE_IAL_pipelines](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) project.


## Relevant Documents
* [OU-SHE Common Testing Environment](https://euclid.roe.ac.uk/projects/sgsshear/wiki/OUSHE_Common_Testing_Environment) on Redmine


## Dependencies

### Internal Euclid Dependencies

* [SHE_PPT 9.6.1](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT)


### External Euclid Dependencies

* [EL_Utils 1.6.2](https://gitlab.euclid-sgs.uk/EuclidLibs/EL_Utils)
* [ST_DataModelTools 9.3.2](https://gitlab.euclid-sgs.uk/ST-DM/ST_DataModelTools)
* [ST_FitsDataModel 9.2.1](https://gitlab.euclid-sgs.uk/ST-DM/ST_FitsDataModel)
* [ST_DataModel 9.2.3](https://gitlab.euclid-sgs.uk/ST-DM/ST_DataModel)
* [Elements 6.2.1](https://gitlab.euclid-sgs.uk/ST-TOOLS/Elements)


### Configuration

**EDEN 3.0**
```
- astropy 5.0
- numpy 1.20.3
- etc
```

### Dependant Projects

* [SHE_IAL_Pipelines](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines)


### Dependant Pipelines

* [Analysis](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) 
* [Calibration](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py)
* [Global Validation](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Global_Validation/PipScript_SHE_Global_Validation.py)
* [Reconciliation](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Reconciliation/PipScript_SHE_Shear_Reconciliation.py)
* [Scaling Experiments](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py)


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
  - [`SHE_CTE_MeasureStatistics`](#she_cte_measurestatistics): For a set of shear measurements, computes statistics on them
  - [`SHE_CTE_MeasureBias`](#she_cte_measurebias): Given a set of statistics files (from `SHE_CTE_MeasureStatistics`), calculates the bias
  - [`SHE_CTE_PlotBias`](#she_cte_plotbias):(Deprecated) Plots the bias as determined by `SHE_CTE_MeasureBias`
  - [`SHE_CTE_PlotPsfSensitivity`](#she_cte_plotpsfsensitivity): (Deprecated)
  - [`SHE_CTE_PrintBias`](#she_cte_printbias): Prints the bias as determined by `SHE_CTE_MeasureBias`
* SHE_CTE_PipelineUtility: Miscellaneous helper executables
  - [`SHE_CTE_DetectorCatalog`](#she_cte_detectorcatalog): Create a sub-catalog of objects within a single detector
  - [`SHE_CTE_ObjectIdSplit`](#she_cte_objectidsplit): Splits a list of objects from the input catalogue into a number of smaller batches
  - [`SHE_CTE_SubObjectIdSplit`](#she_cte_subobjectidsplit): Splits batches of objects into even smaller batches (see `SHE_CTE_ObjectIdSplit`)
  - [`SHE_CTE_ShearEstimatesMerge`](#she_cte_shearestimatesmerge): Merges multiple shear estimate products into one.
  - [`SHE_CTE_CleanupBiasMeasurement`](#she_cte_cleanupbiasmeasurement): Cleans up intermediate files involved in bias measurement
  - [`SHE_CTE_ConvertToHDF5`](#she_cte_converttohdf5): Converts VIS frames and SHE reprojected segmentation maps to HDF5 for more efficient stamp extraction.
  - [`SHE_CTE_UnzipFiles`](#she_cte_unzipfiles): Unzips any gzipped files pointed to by a product, creating a new product pointing to these unzipped files.
  - [`SHE_CTE_ReorderListfiles`](#she_cte_reorderlistfiles): takes a pair of listfiles and reorders the second such that its entries correspond to those in the first.
* SHE_CTE_ScalingExperiments: Executables for a test pipeline to determine how to improve the scaling of the analysis pipeline
  - [`SHE_CTE_SplitFits`](#she_cte_splitfits): splits fits file for an exposure into individual FITS/HDF5 files, one per CCD
  - [`SHE_CTE_CombineSplitFitsListfile`](#she_cte_combinesplitfitslistfile): Combines output from `SHE_CTE_SplitFits` into a single json
  - [`SHE_CTE_ExtractObjects`](#she_cte_extractobjects): For an input catalogue and observation, extracts all the objects from the catalogue that are present in the observation
  - [`SHE_CTE_MakeBatches`](#she_cte_makebatches): Splits the list of objects in the observation into small batches
  - [`SHE_CTE_ExtractStamps`](#she_cte_extractstamps): Given a list of objects, times itself extracting 400x400 pixel postage stamps around the objects
  - [`SHE_CTE_AnalyseRuntime`](#she_cte_analyseruntime): Analyses the runtime of the pipeline, producing graphs
* SHE_CTE_ShearEstimation: Shear Estimation executables
  - [`SHE_CTE_EstimateShear`](#she_cte_estimateshear) *DEPRECATED*: Given a list of objects, estimates their shear
* SHE_CTE_ShearReconcilliation: Shear Reconciliation Executables
  - [`SHE_CTE_ReconcileShear`](#she_cte_reconcileshear): Reconciles different shear measurements for objects from multiple observations.
* SHE_CTE_ShearValidation: Shear validation executables
  - [`SHE_CTE_CrossValidateShear`](#she_cte_crossvalidateshear): (Deprecated)


## Running the software

### `SHE_CTE_MeasureStatistics`

This executable takes the actual and measured shears for each object and from these constructs and stores shear statistics data tables for them (see the LinregressStatistics object in [`SHE_PPT/math.py`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/math.py))

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_MeasureStatistics --details_table <file> --shear_estimates <file> --pipeline_config <file> --she_bias_statistics <file> [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile]
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
| --details_table `<filename>`     | `.xml` intermediate data product pointing to sheSimulatedCatalog FITS tables | yes    | N/A         |
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

**Description:** The filename of an `.xml` data product in the workdir pointing to [sheSimulatedCatalog](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_simulated_catalog.py) FITS tables.

**Source:**  Generated by the `SHE_GST_GenGalaxyImages` executable from [`SHE_GST`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST)


_`shear_estimates`_:

**Description:** The filename of an `.xml` [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html) data product in the workdir. 

**Source:** Generated by `SHE_CTE_EstimateShear`

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
| SHE_CTE_MeasureStatistics_archive_dir `<dir>` | The directory for archiving the data.   | None |
| SHE_CTE_MeasureStatistics_webdav_archive `<True/False>`  | Sets whether to archive to webdav or not | False |
| SHE_CTE_MeasureStatistics_webdav_dir `<dir>` | The webdav directory to archive to | None |

If both these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_bias_statistics`_ 

**Description:** The desired name of the output shear statistics `.xml` data product. This is a sheIntermediateGeneral data product, and points to a number of shear bias statistics [FITS tables](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py), one per estimation method.




_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.








### `SHE_CTE_MeasureBias`

From a set of input bias statistics products, combines these and calculates the bias for each measurement method. Also has a "recovery mode" whereby it will search the workdir for statistics products to use.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_MeasureBias --she_bias_statistics <file>  --pipeline_config <file> --she_bias_measurements <file> [--details_table_head <str>] [--bootstrap-seed <int>] [--recovery_bias_statistics_filename <filename>] [--recovery_bias_measurements_filename <filename>] [--store_measurements_only] [--use_bias_only] [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile]
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

**Description:** The filename of a `.json` listfile pointing to a number of `.xml` sheIntermediateGeneral data products, that point to a shear bias statistics [FITS tables](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py))

**Source:**  Files pointed to in the listfile are generated by the `SHE_CTE_MeasureStatistics` executable (see above). The listfile is generated automatically by the shear calibration pipeline when combining parallel branches of the pipeline together.


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
| SHE_CTE_MeasureBias_archive_dir `<dir>` | The directory for archiving the data.   | None |
| SHE_CTE_MeasureBias_webdav_archive `<True/False>`  | Sets whether to archive to webdav or not | False |
| SHE_CTE_MeasureBias_webdav_dir `<dir>` | The webdav directory to archive to | None |
| SHE_CTE_MeasureBias_number_threads `<int>` | The number of threads (processes) used when reading in the bias statistics files | 8 |


If these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_bias_measurements`_ 

**Description:** The desired name of the output bias measurements `.xml` data product

**Details:** The product is of the type sheIntermediateGeneral, and points to a number of shear bias statistics [FITS tables](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_bias_statistics.py), one per measurement method, which contain the bias measurements in their headers, but also can contain the bias statistics in the body of the table. 


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_PlotBias`

This executable is deprecated, although may be resurrected in the future.

<!-- _**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MeasureStatistics` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_PlotBias --methods <method list>  --root-data-folder <dir> --bias_measurements_head <str> --output_file_name_head <str> [--output_format <str>] [--hide] [--plot_error] [--plot_slopes] [--normed_only] [--unnormed_only] [--workdir <dir>]  [--logdir <dir>] [--profile]
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

Unknown -->








### `SHE_CTE_PlotPsfSensitivity`

This executable is deprecated.

<!-- _**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_PlotPsfSensitivity` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_PlotPsfSensitivity --methods <method list>  --root-data-folder <dir> --bias_measurements_head <str> --output_file_name_head <str> [--output_format <str>] [--hide] [--plot_error] [--normed_only] [--unnormed_only] [--workdir <dir>]  [--logdir <dir>] [--profile]
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

Unknown -->



### `SHE_CTE_PrintBias`

Prints the bias calculated by `SHE_CTE_MeasureBias` out to screen.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_PrintBias` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_PlotPsfSensitivity --she_bias_measurements <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
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

_**Outputs**_
This executable produces no outputs. It only prints values to the screen.


_**Example**_

Unknown









### `SHE_CTE_CleanupBiasMeasurement`

Cleans up intermediate files that are created during the pipeline execution, but are not needed any more.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_CleanupBiasMeasurement` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_CleanupBiasMeasurement --simulation_config <file> --data_images <file> --stacked_data_image <file> --psf_images_and_tables <file> --segmentation_images <file> --stacked_segmentation_image <file> --detections_tables <file> --details_table <file> --shear_estimates <file> --shear_bias_statistics_in <file> --shear_bias_statistics_out <file> --pipeline_config <file>  [--workdir <dir>]  [--logdir <dir>] [--profile]
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

**Description:** A Listfile pointing to detection table ([MER Final Catalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html)) data products

**Source:**  Generated by `SHE_GST_GenGalaxyImages` (See [SHE_GST](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_GST))


_`details_table`_:

**Description:** An `.xml` product pointing to [sheSimulatedCatalog](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_PPT/-/blob/develop/SHE_PPT/python/SHE_PPT/table_formats/she_simulated_catalog.py) FITS tables.

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

This program is designed to be run on intermediate data generated within an execution of the [`SHE Calibration`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data, *_although as its job is to delete files it may throw an exception as the files it has to delete will no longer exist!_* The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.


### `SHE_CTE_DetectorCatalog`

Takes a MER final catalog and the WCS from a single VIS detector and outputs a sub-catalog of only those objects that fall within the footprint of that detector.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_DetectorCatalog` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_DetectorCatalog --workdir <dir> --vis_detector_frame <file> --mer_final_catalog <file> --mer_detector_catalog <file> --she_object_id_list <file> [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**      | **Description**               | **Required** | **Default** |
| :----------------- | :---------------------------- | :----------: | :----------:|
| --workdir `<path>` | Name of the working directory | Yes          | N/A         |
| --logdir `<path>`  | Name of the log directory     | No           | workdir     |

**Input Arguments**

|  **Argument**                       | **Description**                                         | **Required** | **Default** |
| :---------------------------------- | :------------------------------------------------------ | :----------: | :----------:|
| --vis_detector_frame `<filename>`   | visCalibratedFrame data product containing one detector | Yes          | N/A         |
| --mer_final_catalog `<filename>`    | merFinalCatalog data product                            | Yes          | N/A         |

**Output Arguments**

|  **Argument**                       | **Description**                                                      | **Required** | **Default** |
| :---------------------------------- | :------------------------------------------------------------------- | :----------: | :----------:|
| --mer_detector_catalog `<filename>` | merFinalCatalog containing objects within the VIS detector footprint | Yes          | N/A         |
| --she_object_id_list `<filename>`   | sheObjectIdList containing the objects IDs of those objects          | Yes          | N/A         |

_**Inputs**_

_`vis_detector_frame`_:

**Description:** The filename of a VIS calibrated fram data product that has been preprocessed to contain data for only a single detector.

**Source:** SHE_CTE_DetectorSplit

_`mer_final_catalog`_:

**Description:** The filename of a merFinalCatalog data product.

**Source:** MER

_**Outputs**_

_`mer_detector_catalog`_

**Description** A merFinalCatalog data product containing only those object from the input catalog that fall within the footprint of the input VIS detector.

_`she_object_id_list`_

**Description** A sheObjectIdList data product containing the object IDs of objects from the input MER catalog that fall within the footprint of the input VIS detector.


### `SHE_CTE_ObjectIdSplit`

Finds objects in the MER catalogue that are in the observation, then divides them into batches for processing by other executables. This is part of the analysis pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ObjectIdSplit` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ObjectIdSplit --mer_final_catalog_tables <file> --data_images <file>  --pipeline_config <file> --object_ids <file> --batch_mer_catalogs <file> [--skip_vis_non_detections] [--workdir <dir>]  [--logdir <dir>] [--profile]
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
| --skip_vis_non_detections  | If set, will not include objects with VIS_DET==0 in the output catalogues  | no           | false       |


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
| SHE_CTE_ObjectIdSplit_batch_ids `<str>` | A list of object IDs to use | None |
| SHE_CTE_ObjectIdSplit_grouping_radius `<float>` | The grouping radius (in arcseconds) | 1.0 |


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

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.





### `SHE_CTE_SubObjectIdSplit`

This executable is virtually identical to SHE_CTE_ObjectIdSplit, save for its pipeline_config options. It further splits a batch of objects (from SHE_CTE_ObjctIDSplit) into smaller sub-batches.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_SubObjectIdSplit` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_SubObjectIdSplit --mer_final_catalog_tables <file> --data_images <file>  --pipeline_config <file> --object_ids <file> --batch_mer_catalogs <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
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
| --batch_mer_catalogs `<filename>`  | Name of the `.json` listfile containing the list of merFinalCatalog products for each sub-batch | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`mer_final_catalog_tables`_:

**Description:** The filename of a `.json` listfile pointing to an `.xml` [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) data product from a batch.

**Source:** This is produced by SHE_CTE_ObjectIdSplit


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
| SHE_CTE_SubObjectIdSplit_batch_ids `<str>` | A list of object IDs to use | None |
| SHE_CTE_ObjectIdSplit_grouping_radius `<float>` | The grouping radius (in arcseconds) *note: this is the same option as for SHE_CTE_ObjectIdSplit* | 1.0 |


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`object_ids`_ 

**Description** The desired name for a listfile pointing to numerous [sheObjectIdList](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_objectidlist.html) products, one per sub-batch.

_`batch_mer_catalogs`_ 

**Description** The desired name for a listfile pointing to numerous [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) products, one per batch. These contain the same information as the input MER final catalogues, just only the rows for each object in the sub-batch.

_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_ShearEstimatesMerge`

Given input lists shear estimates and LensMC chains, combines them into two tables, one for the estimates, and one for the chains. This is a merge point executable in the Analysis pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ShearEstimatesMerge` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ShearEstimatesMerge --shear_estimates_product_listfile <file> --she_lensmc_chains_listfile <file>  --pipeline_config <file> ----merged_she_measurements <file> --merged_she_lensmc_chains <file> [--workdir <dir>]  [--logdir <dir>] [--profile]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE Analysis`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.



### `SHE_CTE_ConvertToHDF5`

Given input VIS Calibrated frames and SHE reprojected segmentation maps, converts these to a HDF5 file for more efficient stamp extraction.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ConvertToHDF5` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ConvertToHDF5 --vis_frame_prod <file> --remapped_seg_prod <file>  --output_hdf5 <file> [--workdir <dir>]  [--logdir <dir>]
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
| --vis_frame_prod `<filename>`     | DpdVisCalibratedFrame product | yes    | N/A         |
| --remapped_seg_prod`<filename>`     | DpdSheExposureReprojectedSegmentationMap product | yes    | N/A         |



**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_hdf5 `<filename>`  | Name of the `.h5` file | yes          | N/A |


**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --chunksize     | The size of chunk (in pixels) for the data  | no           | 100       |


_**Inputs**_

_`vis_frame_prod`_:

**Description:** The filename [visCalibratedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_calibratedframe.html) data product.

**Source:** This is provided as input to the Analysis pipeline


_`remapped_seg_prod`_:

**Description:** The filename of a [sheExposureReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_exposurereprojectedsegmentationmap.html) data product.

**Source:** This is produced by SHE_MER_RemapMosaic 



_**Outputs**_

_`output_hdf5`_ 

**Description** The desired name for the output HDF5 file.



### `SHE_CTE_UnzipFiles`

Given an input data product, unzipps all gzipped files it points to, and creates a new product pointing to these unzipped files. Any files that are not zipped will be left alone.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_UnzipFiles` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_UnzipFiles --input_prod <file> --output_prod <file> [--workdir <dir>]  [--logdir <dir>]
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
| --input_product `<filename>`     | Any xml data product | yes    | N/A         | 


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --output_product `<filename>`  | Identical product to the input, but pointing to unzipped files | yes          | N/A |


**Options**

None


_**Inputs**_

_`input_product`_:

**Description:** The filename of any xml data product pointing to files which are gziped.
**Source:** Any


_**Outputs**_

_`output_product`_ 

**Description** The desired name for the output product. This will be the same product type as the input product, and contain the same contents, save for the files it links to will be the unzipped versions of these files.

### `SHE_CTE_ReorderListfiles`

Given a reference listfile and an input one, will reorder the input listfile such that its entries correspond to those in the reference listfile, according to user supplied matching criteria.

For example, say we have as a reference listfile containing a list of products with a quantity, where the quantities for each product in the list are [a, b, c, d]. Then consider we have an input listfile with products containing the same quantity, but in different order, e.g. [d, c, b, a]. This executable will re-order the input listfile to have entries in the order such that the quantity appears [a, b, c, d]. E,g, the nth item in the reference listfile corresponds to the nth item in the output listfile.

This executable can also be used to expand a listfile. So, for example, if we have a reference list of exposures (belonging to several observations, but each observation may be represented more than once) such that some observations are duplicated, and a list of input products where there is only one per observation, then the executable can be used to expand the input listfile so there is one entry per reference exposure. E.g., if the reference exposures have observations [a, a, b, b] and the input listfile has observations [a, b] then the output listfile would have observations [a, a, b, b].

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ReorderListfiles` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ReorderListfiles --reference_listfile <file> --input_listfile <file> --output_listfile <file> --reference_path <XML path> --input_path <XML path> [--allow_duplicates] [--workdir <dir>]  [--logdir <dir>]
```
with the following options:

**Common Elements Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --workdir `<path>`          | Name of the working directory, where input data is stored and output data will be created. | no          | .         |
| --logdir `<path>`     | Name of the directory for misc log files, e.g. profiling files. | no| . |


**Input Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --reference_listfile `<filename>`     | Any listfile pointing to products| yes    | N/A         | 
| --input_listfile `<filename>`     | Any listfile pointing to products| yes    | N/A         | 


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --reference_listfile `<filename>`     | Any listfile pointing to products| yes    | N/A         | 


**Options**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --allow_duplicates     | If set, will allow the reference listfile's products to contain duplicate values, and the output listfile will be expended to include duplicate items where necessary | No   | Unset        | 


_**Inputs**_

_`reference_listfile`_:

**Description:** The filename of a listfile pointing to products which the user intends to use for reference. The input listfile will be reordered according to the entries in the reference listfile
**Source:** Any

_`input_listfile`_:

**Description:** The filename of a listfile pointing to products which the user intends to have reordered relative to the reference listfile
**Source** Any

_`reference_path`_:

**Description:** The xml path within the reference product to the quantity that the user wishes to reorder relative to. E.g. `Data.ObservationSequence.ObservationId`

_`input_path`_:

**Description:** The xml path within the input product to the quantity that the user wishes to reorder relative to. E.g. `Data.ObservationSequence.ObservationId`

_**Options**_:

_`allow_duplicates`_:

**Description** If set, then the executable will allow the reference listfile to contain duplicate values, and it may produce an output listfile with duplicate entries.

_**Outputs**_

_`output_listfile`_ 

**Description** The desired name for the output listfile, which is the reordered version of the input listfile


### `SHE_CTE_SplitFits`

Takes the input exposure and splits it into individual files, one per CCD. This executable can produce HDF5 or FITS files depending upon the configuration options set. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_SplitFits` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_SplitFits --input_fits_json <file> --pipeline_config <file> --output_json <file> --timing_info <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.




### `SHE_CTE_CombineSplitFitsListfile`

Combines the output from SHE_CTE_SplitFits into a single json. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_CombineSplitFitsListfile` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_CombineSplitFitsListfile --input_listfile <file> --pipeline_config <file> --output_json <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.



### `SHE_CTE_ExtractObjects`

Given an input catalogue of objects and an observation, extracts the objects from the catalogue that are in the observation, recording which CCD the object is on for each exposure. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ExtractObjects` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ExtractObjects --stacked_image <file> --exposures <file> --catalogue_listfile <file> --pipeline_config <file> --output_objects_list <file> --combined_catalogue <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.




### `SHE_CTE_MakeBatches`

Takes the list of all the objects in the observation and splits them into small batches. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_MakeBatches` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_MakeBatches --objects_list <file> --pipeline_config <file> --batch_listfile <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.







### `SHE_CTE_ExtractStamps`

Given the input list of objects, extracts 400x400 pixel postage stamps from each CCD the object is in and times how long it takes to do this. This is part of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ExtractStamps` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ExtractStamps --batch_file <file> --exposures <file> --pipeline_config <file> --timing_info <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.







### `SHE_CTE_AnalyseRuntime`

Analyses the runtime of SHE_CTE_SplitFITS and SHE_CTE_ExtractStamps and produces a number of graphs and a json file with some statistics. This is the final step of the ScalingExperiments pipeline, used to determine the optimal method for reading postage stamps in at scale.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_AnalyseRuntime` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_AnalyseRuntime --stamp_timing_listfile <file> --split_timing_listfile <file> --pipeline_config <file> --results <file> [--workdir <dir>]  [--logdir <dir>]
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

This program is designed to be run on intermediate data generated within an execution of the [`SHE ScalingExperiments`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Scaling_Experiments/PipScript_SHE_Scaling_Experiments.py) pipeline. Please see the documentation of that pipeline for an example run. After that pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.








### `SHE_CTE_EstimateShear` *DEPRECATED*

This executable measures the shear for a list of input objects. It is used in the analysis pipeline, as well as the calibration pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_EstimateShear` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_EstimateShear --data_images <file> --stacked_image <file> --psf_images_and_tables <file> --segmentation_images <file> stacked_segmentation_image <file> --detections_tables <file> --object_ids <file> --ksb_training_data <file> --lensmc_training_data <file> -momentsml_training_data <file> --regauss_training_data <file> --mdb <file> --pipeline_config <file>  --shear_estimates_product <file> --she_lensmc_chains <file>  [--methods <str>] [--chains_method <str>] [--memmap_images] [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile] [--debug]
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
| --data_images `<filename>`     | `.json` listfile pointing to a series of data image products | yes    | N/A         |
| --stacked_image `<filename>`     | `.xml` data product for the stacked image | yes    | N/A         |
| --psf_images_and_tables `<filename>`     | `.json` listfile pointing to psf image products | yes    | N/A         |
| --segmentation_images `<filename>`     | `.json` listfile pointing to a series of segmentation map data products | yes    | N/A         |
| --stacked_segmentation_image `<filename>`     | `.xml` data product for the segmentation map for the stacked image | yes    | N/A         |
| --detections_tables `<filename>`     | `.json` listfile pointing to a series of detections table products | yes    | N/A         |
| --object_ids `<filename>`     | `.xml` data product containing the list of objects to have their shear measured | yes    | N/A         |
| --ksb_training_data `<filename>`     | `.xml` data product containing KSB training data | yes    | N/A         |
| --lensmc_training_data `<filename>`     | `.xml` data product containing LensMC training data | yes    | N/A         |
| --momentsml_training_data `<filename>`     | `.xml` data product containing MomentsML training data | yes    | N/A         |
| --regauss_training_data `<filename>`     | `.xml` data product containing REGAUSS training data | yes    | N/A         |
| --mdb `<filename>`     | `.xml` mission database data product | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --shear_estimates_product `<filename>`  | `.xml` data product for the shear estimates | yes          | N/A |
| --she_lensmc_chains `<filename>`  | `.xml` data product for the LensMC chains | yes          | N/A |


**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |
| --debug (`store true`)    | If set, will only process the first 1000 objects in the input object list  | no           | false       |
| --methods `<str>`    | The list of shear estimation methods to use. If blank, will use all the methods. The methods to choose from are "KSB", "REGAUSS", "LensMC" and "MomentsML"  | no           | None       |
| --chains_method `<str>`    | Which shear estimation method to generate chains with. The methods to choose from are "KSB", "REGAUSS", "LensMC" and "MomentsML"  | no           | "LensMC"      |
| --memmap_images (`store true`)    | if set, opens the images with memmap | no           | false       |
| --lensmc_stamp_size `<int>`    | The requested stamp size in pixels | no           | 384       |
| --lensmc_x_buffer `<int>`    | Do not fit object if closer to the edge of the detector than this number of pixels (x direction) | no           | 3       |
| --lensmc_y_buffer `<int>`    | Do not fit object if closer to the edge of the detector than this number of pixels (y direction) | no           | 3       |
| --lensmc_no_mask_dilation `<bool>`    | do not dilate mask by 1 pixel | no           | false       |
| --lensmc_hl_to_exp `<float>`    |Half-light radius of the bulge to exponential scalelength of the disc | no           | 0.15       |
| --lensmc_n_bulge `<float>`    | Bulge Sersic index; available: n=(1., 1.5, 2., 2.5, 3., 3.5, 4.) | no           | 1       |
| --lensmc_n_disc `<float>`    | Disc Sersic index; available: n=1 | no           | 1       |
| --lensmc_e_max `<float>`    | Hard upper bound on ellipticity | no           | 0.99       |
| --lensmc_re_max `<float>`    | Hard upper bound on effective radius | no           | 2.0       |
| --lensmc_delta_max `<float>`    | Hard upper bound on position offset | no           | 0.3       |
| --lensmc_e_flag `<float>`    | Flagging threshold for ellipticity | no           | 0.98       |
| --lensmc_re_flag `<float>`    | Flagging threshold for effective radius | no           | 1.8       |
| --lensmc_delta_flag`<float>`    | Flagging threshold for position offset | no           | 0.28       |
| --lensmc_disc_only `<bool>`    | Whether to fit only for a disc component | no           | False       |
| --lensmc_psf_oversampling `<int>`    | PSF oversampling factor | no           | 5       |
| --lensmc_seed `<int>`    | Seed the random sampler (-1 means use the has of the object's ID) | no           | -1       |
| --lensmc_shape_noise `<float>`    | Shape noise standard deviation if not provided by training data | no           | 0.25       |
| --lensmc_return_chains `<bool>`    | Whether to return the chains | no           | false       |
| --lensmc_fast_mode `<bool>`    | Enable LensMC fast mode. Override any sampling settings and produce MAP estimate (without MCMC/errors/intcal) | no           | false       |
| --lensmc_include_vis_undetected `<bool>`    | Measure all objects, including those that have not been detected by VIS | no           | false       |
| --lensmc_monitor `<bool>`    | Monitor data and visually check consistency of input data | no           | false       |


_**Inputs**_

_`data_images`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [visCalibratedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_calibratedframe.html) data products.

**Source:**  This is provied as input to the analysis pipeline


_`stacked_image`_:

**Description** The filename of a [visStackedFrame](https://euclid.esac.esa.int/dm/dpdd/latest/visdpd/dpcards/vis_visstackedframe.html) data product.

**Source** This is provied as input to the analysis pipeline


_`psf_images_and_tables`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [shePSFModelImage](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_psfmodelimage.html) data products.

**Source:**  This is produced by SHE_PSFToolkit_ModelPSFs


_`segmentation_images`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [sheExposureReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_exposurereprojectedsegmentationmap.html) data products.

**Source:**  This is produced by SHE_MER_RemapMosaic


_`stacked_segmentation_image`_:

**Description:** The filename of a [sheStackReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_stackreprojectedsegmentationmap.html) data product.

**Source:**  This is produced by SHE_MER_RemapMosaic


_`detections_tables`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) data products.

**Source:**  This can be produced in a number of different ways depending upon where this executable is called, but in the simplest case, it is provided as an input to the pipeline. Another source is from SHE_CTE_(Sub)ObjectIdSplit, which takes the pipeline's input detections tables and extracts only the relevant rows from them.


_`object_ids`_:

**Description:** The filename of a [sheObjectIdList](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_objectidlist.html) data product.

**Source:**  This is produced by SHE_CTE_SubObjectIdSplit in the analysis pipeline.


_`ksb_training_data`_:

**Description:** The filename of a [sheKsbTraining](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_ksbtraining.html) data product.

**Source:**  This is provided as input to the pipeline


_`lensmc_training_data`_:

**Description:** The filename of a [sheLensMcTraining](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmctraining.html) data product.

**Source:**  This is provided as input to the pipeline


_`momentsml_training_data`_:

**Description:** The filename of a [sheMomentsMlTraining](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_momentsmltraining.html) data product.

**Source:**  This is provided as input to the pipeline


_`regauss_training_data`_:

**Description:** The filename of a [sheRegaussTraining](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_regausstraining.html) data product.

**Source:**  This is provided as input to the pipeline


_`mdb`_:

**Description:** The filename of the `.xml` [mission database](https://euclid.roe.ac.uk/projects/missiondatabase/wiki/Wiki).

**Source:**  This is provided as input to the pipeline


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default values** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_EstimateShear_methods `<str>` | A list of shear estimation methods to use   | `["KSB","REGAUSS","MomentsML","LensMC"]` |
| SHE_CTE_EstimateShear_chains_method `<str>`  | The method to generate chains for. The methods to choose from are "KSB", "REGAUSS", "LensMC" and "MomentsML". | `"LensMC"` |
| SHE_CTE_EstimateShear_memmap_images `<True/False>` | Use memmap for opening images or not| False |
| SHE_LensMC_stamp_size `<int>`    | The requested stamp size in pixels  | 384       |
| SHE_LensMC_x_buffer `<int>`    | Do not fit object if closer to the edge of the detector than this number of pixels (x direction)           | 3       |
|  SHE_LensMC_y_buffer `<int>`    | Do not fit object if closer to the edge of the detector than this number of pixels (y direction)          | 3       |
|  SHE_LensMC_no_mask_dilation `<bool>`    | do not dilate mask by 1 pixel          | false       |
|  SHE_LensMC_hl_to_exp `<float>`    |Half-light radius of the bulge to exponential scalelength of the disc          | 0.15       |
|  SHE_LensMC_n_bulge `<float>`    | Bulge Sersic index; available: n=(1., 1.5, 2., 2.5, 3., 3.5, 4.)           | 1       |
|  SHE_LensMC_n_disc `<float>`    | Disc Sersic index; available: n=1          | 1       |
|  SHE_LensMC_e_max `<float>`    | Hard upper bound on ellipticity         | 0.99       |
|  SHE_LensMC_re_max `<float>`    | Hard upper bound on effective radius           | 2.0       |
|  SHE_LensMC_delta_max `<float>`    | Hard upper bound on position offset         | 0.3       |
|  SHE_LensMC_e_flag `<float>`    | Flagging threshold for ellipticity       | 0.98       |
|  SHE_LensMC_re_flag `<float>`    | Flagging threshold for effective radius        | 1.8       |
|  SHE_LensMC_delta_flag`<float>`    | Flagging threshold for position offset       | 0.28       |
|  SHE_LensMC_disc_only `<bool>`    | Whether to fit only for a disc component       | False       |
|  SHE_LensMC_oversampling `<int>`    | PSF oversampling factor      | 5       |
|  SHE_LensMC_seed `<int>`    | Seed the random sampler (-1 means use the hash of the object's ID)       | -1       |
|  SHE_LensMC_shape_noise `<float>`    | Shape noise standard deviation if not provided by training data       | 0.25       |
|  SHE_LensMC_return_chains `<bool>`    | Whether to return the chains        | false       |
|  SHE_LensMC_fast_mode `<bool>`    | Enable LensMC fast mode. Override any sampling settings and produce MAP estimate (without MCMC/errors/intcal)     | false       |
|  SHE_LensMC_include_vis_undetected `<bool>`    | Measure all objects, including those that have not been detected by VIS       | false       |
|  SHE_LensMC_monitor `<bool>`    | Monitor data and visually check consistency of input data        | false       |

If both these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`shear_estimates_product`_ 

**Description** The desired name of the output shear measurements `.xml` data product

**Details** The product is of the type [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html).


_`she_lensmc_chains`_ 

**Description** The desired name of the output LensMC chains `.xml` data product

**Details** The product is of the type [sheLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmcchains.html).


_**Example**_

This program is designed to be run on intermediate data generated within an execution of the [Calibration](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Calibration/PipScript_SHE_Shear_Calibration.py) and [Analysis](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Analysis/PipScript_SHE_Shear_Analysis.py) pipelines. Please see the documentation for these pipelines for an example run. After the pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_ReconcileShear`

This executable takes all the shear measurements for objects (from different observations) for a given MER tile, and reconciles any different measurements for each individual object. This is the main executable in the Reconciliation pipeline.

_**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_ReconcileShear` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_ReconcileShear --she_validated_measurements_listfile <file> --she_lensmc_chains_listfile <file> --mer-final-catalog <file> --she_reconciliation_config <file> --pipeline_config <file> --she_reconciled_measurements <file> --she_reconciled_lensmc_chains <file>  [--method <str>] [--chains_method <str>] [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile] [--debug]
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
| --she_validated_measurements_listfile `<filename>`     | `.json` listfile pointing to a series validated measurements | yes    | N/A         |
| --she_lensmc_chains_listfile `<filename>`     | `.json` listfile pointing to a series of LensMC chains | yes    | N/A         |
| --psf_images_and_tables `<filename>`     | `.json` listfile pointing to psf image products | yes    | N/A         |
| --mer_final_catalog `<filename>`     | `.xml` data product for a MER tile | yes    | N/A         |
| --she_reconciliation_config `<filename>`     | `.xml` data product for the reconciliation configuration | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --she_reconciled_measurements `<filename>`  | `.xml` data product for the reconciled measurements | yes          | N/A |
| --she_reconciled_lensmc_chains `<filename>`  | `.xml` data product for reconciled LensMC chains | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --method `<str>`   | Which reconciliation method to use. The options are "InvVar", "Best" and "ShapeWeight". | no           | `"InvVar"`       |
| --chains_method `<str>`   | Which LensMc chains reconciliation method to use. The options are "InvVar", "Best", "ShapeWeight" and "keep" | no           | `"keep"`     |
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`she_validated_measurements_listfile`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [sheValidatedMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_validatedmeasurements.html) data products.

**Source:**  This is provied as an input to the pipeline


_`she_lensmc_chains_listfile`_:

**Description** The filename of a `.json` Listfile pointing to a number of [sheLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmcchains.html) data products.

**Source** This is the output from SHE_CTE_EstimateShear


_`mer_final_catalog`_:

**Description:** The filename of a [merFinalCatalog](https://euclid.esac.esa.int/dm/dpdd/latest/merdpd/dpcards/mer_finalcatalog.html) data product.

**Source:**  This is provided as an input to the pipeline


_`segmentation_images`_:

**Description:** The filename of a `.json` Listfile pointing to a number of [sheExposureReprojectedSegmentationMap](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_exposurereprojectedsegmentationmap.html) data products.

**Source:**  This is produced by SHE_MER_RemapMosaic


_`she_reconciliation_config`_:

**Description:** The filename of a sheReconciliationConfig (no link available) data product.

**Source:**  This is provided as input to the pipeline


_`pipeline_config`_:

**Description:** One of the following:

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.

The `.txt` pipeline configuration file may have any number of configuration arguments which apply to other executables, in addition to optionally any of the following which apply to this executable:


|  **Options**                | **Description**                                                       | **Default values** |
| :------------------------   | :-------------------------------------------------------------------- | :----------:|
| SHE_CTE_ReconcileMeasurements_method `<str>` | The reconciliation method to use. The options are "InvVar", "Best" and "ShapeWeight"  | `"InvVar"` |
| SHE_CTE_ReconcileMeasurements_chains_method `<str>`  | The LensMC reconciliation method to use. The options are "InvVar", "Best", "ShapeWeight" and "keep"  | `"keep"` |

If both these arguments are supplied in the pipeline configuration file and the equivalent command-line arguments are set, the command-line arguments will take precedence.

**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`she_reconciled_measurements`_ 

**Description** The desired name of the output shear reconciled measurements `.xml` data product

**Details** The product is of the type [sheReconciledMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_reconciledmeasurements.html).


_`she_reconciled_lensmc_chains`_ 

**Description** The desired name of the output LensMC chains `.xml` data product

**Details** The product is of the type [sheReconciledLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_reconciledlensmcchains.html).


_**Example**_

This executable is part of the [Reconciliation pipeline](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines/-/blob/develop/SHE_Pipeline/auxdir/SHE_Shear_Reconciliation/PipScript_SHE_Shear_Reconciliation.py). Please see the documentation for that pipeline for an example run. After the pipeline has been run once, this program can be re-run on the generated intermediate data. The command used for the execution of this program will be stored near the top of the log file for its original execution, which can be found in logdir within the workdir after execution.






### `SHE_CTE_CrossValidateShear`

This executable is deprecated.

<!-- _**Running the Program on EDEN/LODEEN**_

To run `SHE_CTE_CrossValidateShear` on Elements use the following command:
```bash
E-Run SHE_CTE 9.3 SHE_CTE_CrossValidateShear --shear_estimates_product <file> --shear_lensmc_chains <file> --pipeline_config <file> --cross_validated_shear_estimates_product <file> [--primary_method <str>] [--workdir <dir>]  [--logdir <dir>] [--archive_dir <dir>] [--webdav_dir <dir>] [--webdav_archive] [--profile] [--dry_run]
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
| --she_estimates_product `<filename>`     | `.xml` shear estimates data product | yes    | N/A         |
| --she_lensmc_chains `<filename>`     | `.xml` LensMC chains data product | yes    | N/A         |
| --pipeline_config `<filename>` | Config file for the pipeline, containing a number of key-value pairs | no          | None         |


**Output Arguments**

|  **Argument**               | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --cross_validated_shear_estimates_product `<filename>`  | `.xml` data product for the cross validated shear measurements | yes          | N/A |

**Options**

|  **Options**                | **Description**                                                       | **Required** | **Default** |
| :------------------------   | :-------------------------------------------------------------------- | :----------: | :----------:|
| --primary_method `<str>`   | Which estimation method to consider the primary one, and compare against the other methods | no           | `"LensMC"`       |
| --dry_run (`store true`)   | If set, will not actually process any data | no           | `False`     |
| --profile (`store true`)    | If set, will generate profiling information via cProfile in the logdir  | no           | false       |


_**Inputs**_

_`she_estimates_product`_:

**Description:** The filename of an `.xml` [sheMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_measurements.html) data product.

**Source:**  This is produced by SHE_CTE_EstimateShear


_`she_lensmc_chains`_:

**Description** The filename of an `.xml` [sheLensMcChains](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_lensmcchains.html) data product.

**Source** This is the output from SHE_CTE_EstimateShear


_`pipeline_config`_:

**Description:** One of the following:

This executable takes `pipeline_config` as an input, however it does not use any of the config options contained within the file.

1. The word "None" (without quotes), which signals that default values for all configuration parameters shall be used.
1. The filename of an empty `.json` listfile, which similarly indicates the use of all default values.
1. The filename of a `.txt` file in the workdir listing configuration parameters and values for executables in the current pipeline run. This shall have the one or more lines, each with the format "SHE_MyProject_config_parameter = config_value".
1. The filename of a `.xml` data product of format DpdSheAnalysisConfig, pointing to a text file as described above. The format of this data product is described in detail in the Euclid DPDD at https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_calibrationconfig.html.
1. The filename of a `.json` listfile which contains the filename of a `.xml` data product as described above.

Any of the latter three options may be used for equivalent functionality.


**Source:** One of the following:

1. May be generated manually, creating the `.txt` file with your text editor of choice.
1. Retrieved from the EAS, querying for a desired product of type DpdSheAnalysisConfig.
1. If run as part of a pipeline triggered by the [`SHE_Pipeline_Run`](https://gitlab.euclid-sgs.uk/PF-SHE/SHE_IAL_Pipelines) helper script, may be created automatically by providing the argument `--config_args ...` to it (see documentation of that executable for further information).

_**Outputs**_

_`cross_validated_shear_estimates_product`_ 

**Description** The desired name of the output cross validated shear measurements `.xml` data product

**Details** The product is of the type [sheValidateddMeasurements](https://euclid.esac.esa.int/dm/dpdd/latest/shedpd/dpcards/she_validatedmeasurements.html).


_**Example**_

?? Output from whcih pipeline? -->






## Troubleshooting


### A test failed when I ran "make test"

_**Ensure you have the most up-to-date version of the project and all its dependencies**_

It's possible the issue you're hitting is a bug that's already been fixed, or could be due to locally-installed versions of projects on the develop branch no longer being compatible with a newly-deployed version of another dependency on CODEEN. If you're running on the develop branch and have installed locally, pull the project, call `make purge`, and install again, and repeat for all dependencies you've installed locally. Try running `make test` again to see if it works.

_**Report the failing test to the developers**_

If the test still fails, please report it to the active developers listed above, ideally by creating a GitLab issue, or else by e-mailing them.

_**Try running the desired code**_

Tests can fail for many reasons, and a common reason is that the code is updated but not the test. This could lead to the test failing even if the code works properly. After you've reported the issue, you can try to run the desired code before the issue with the failing test has been fixed. There's a decent chance that the bug might only be in the test code, and the executable code will still function.

### An exception was raised, what do I do?

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
