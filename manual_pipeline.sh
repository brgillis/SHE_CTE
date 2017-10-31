#!/usr/bin/bash

WORKDIR="/home/user/euclid-ial/workspace/tests/she_analysis"
LOGDIR=$WORKDIR"/logdir"

DATA_IM="sim_data_images.json"
PSF_CAL="mock_psf_calibration_products.json"
SEG_IM="mock_segmentation_products.json"
DAL_TAB="mock_details_tables.json"
DTC_TAB="mock_detections_tables.json"
AST_PROD="mock_astrometry_products.json"
AOCS_PROD="mock_aocs_time_series_products.json"
MT_PROD="mock_mission_time_products.json"
GP_TAB="mock_galaxy_population_priors_table.fits"
CAL_PROD="mock_calibration_parameters_product.bin"
CAL_LF="mock_calibration_parameters_listfile.json"
SV_TAB="mock_shear_validation_statistics_table.fits"
PSF_IMTAB="mock_psf_images_and_tables.fits"
SE_PROD="mock_shear_estimates.bin"
SE_LF="mock_shear_estimates.json"
VSE_TAB="mock_validated_shear_estimates.fits"

mkdir -p $WORKDIR
mkdir -p $LOGDIR

CMD="E-Run SHE_CTE 0.2 SHE_CTE_MakeMockAnalysisData --data_images $DATA_IM --psf_calibration_products $PSF_CAL --segmentation_images $SEG_IM --detections_tables $DTC_TAB --astrometry_products $AST_PROD --aocs_time_series_products $AOCS_PROD --mission_time_products $MT_PROD --galaxy_population_priors_table $GP_TAB --calibration_parameters_product $CAL_PROD --calibration_parameters_listfile $CAL_LF --shear_validation_statistics_table $SV_TAB --workdir $WORKDIR --logdir $LOGDIR"

echo $CMD
`$CMD`

if [ $? -ne 0 ]; then
    exit
fi

CMD="E-Run SHE_GST 1.2 GenGalaxyImages --config-file /home/user/Work/Projects/SHE_GST/SHE_GST_GalaxyImageGeneration/conf/SHE_GST_GalaxyImageGeneration/SampleStamps.conf --dithering_scheme 2x2 --data_images $DATA_IM --segmentation_images $SEG_IM --detections_tables $DTC_TAB --details_tables $DAL_TAB  --psf_images_and_tables $PSF_IMTAB --workdir $WORKDIR --logdir $LOGDIR"

echo $CMD
`$CMD`

if [ $? -ne 0 ]; then
    exit
fi

CMD="E-Run SHE_CTE 0.2 SHE_CTE_FitPSFs --data_images $DATA_IM --detections_tables $DTC_TAB --astrometry_products $AST_PROD --aocs_time_series_products $AOCS_PROD --mission_time_products $MT_PROD --psf_calibration_products $PSF_CAL --psf_images_and_tables $PSF_IMTAB --workdir $WORKDIR --logdir $LOGDIR"

echo $CMD
`$CMD`

if [ $? -ne 0 ]; then
    exit
fi

CMD="E-Run SHE_CTE 0.2 SHE_CTE_EstimateShear --data_images $DATA_IM --psf_images_and_tables $PSF_IMTAB --segmentation_images $SEG_IM --detections_tables $DTC_TAB --galaxy_population_priors_table $GP_TAB --calibration_parameters_product $CAL_PROD --calibration_parameters_listfile $CAL_LF --shear_estimates_product $SE_PROD --shear_estimates_listfile $SE_LF --workdir $WORKDIR --logdir $LOGDIR"

echo $CMD
`$CMD`

if [ $? -ne 0 ]; then
    exit
fi

CMD="E-Run SHE_CTE 0.2 SHE_CTE_ValidateShear --shear_estimates_product $SE_PROD --shear_estimates_listfile $SE_LF --shear_validation_statistics_table $SV_TAB --validated_shear_estimates_table $VSE_TAB --workdir $WORKDIR --logdir $LOGDIR"

echo $CMD
`$CMD`

