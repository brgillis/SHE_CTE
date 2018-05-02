#!/bin/bash

./build.x86_64-co7-gcc48-o2g/run SHE_CTE_EstimateShear --data_images CalibratedFrames.json --stacked_image DpdVisStackedFrame.xml --psf_images_and_tables EUC_SHE_PSF_IMTAB.json --segmentation_images MosaicFrames.json --stacked_segmentation_image EUC_SHE_SHE_Mosaic_Stack_2018315T15328.0Z_00.00.fits --detections_tables DetectionsCatalogs.json --shear_estimates_product EUC_SHE_SHEAR_ESTIMATES.xml --workdir /home/brg/Data/sc3-workdir --logdir /home/brg/Data/sc3-workdir/logdir --debug
