#!/bin/bash

./build.x86_64-co7-gcc48-o2g/run SHE_CTE_FitPSFs --data_images CalibratedFrames.json --detections_tables DetectionsCatalogs.json --psf_images_and_tables EUC_SHE_PSF_IMTAB.json --workdir /home/brg/Data/sc3-workdir --logdir /home/brg/Data/sc3-workdir/logdir
