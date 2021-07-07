# NICUWatch Testing

## Introduction

This suite is designed to test NICU scripts and workflow processing in order to identify issues/bugs which may arise due to development or connection issues and limit them during true project processing.

The suite is run by using a nextflow test script designed to run all `NICUWatch` steps.

Default project `T101-101-NEO-Cyberdyne` will be used for all testing.

## Required Files

These are the current files required to run the test suite, which must be run from the `~/NICUWatch/test` directory.

```
# Required files and structure.

$> tree ~/NICUWatch/test
.
├── manifest
│   └── NeoSeq_Test_Manifest.xlsx
├── neoseq.test.nf
└── README.md

$> tree /scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00000000/Staging
├── 20329101160_S1_R1_001.fastq.gz.md5
├── 20329101160_S1_R2_001.fastq.gz.md5
├── 20329120647_S2_R1_001.fastq.gz.md5
├── 20329120647_S2_R2_001.fastq.gz.md5
├── 20329120769_S3_R1_001.fastq.gz.md5
└── 20329120769_S3_R2_001.fastq.gz.md5
```

## Running Test Suite

Step by step instructions on how to run a test.

```
# Move to the NICUWatch test directory.
$> cd ~/NICUWatch/test

# Create a screen session (suggested)
$> screen -S neoseq-test

# Load needed modules
$> ml ucgd_modules/dev

# Launch the test workflow
$> nextflow run neoseq.test.nf

```

## Test Run Overview.

* The workflow starts by uploading `NeoSeq_Test_Manifest.xlsx` to the `UCGD_Team/NeoSeq/NeoSeq_manifest_final` directory.

* The `NeoSeq_Test_Manifest.xlsx` file is picked up by NICUWatch and ran as if it was a new project.

NICUWatch is run in `--test` mode for the above steps which will use the fabric dev environment and funnel all SNS messages to a single point/user.


## Important post-processing steps.

After the test suite is completed a few steps are necessary to reset for future use.  

* Run the following command still in the screen session:

```
rm -rf work/ .nextflow*
```

* Move the test fastq from the test directory back to the staging.

```
$> mv /scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00000000/T101-101-NEO-Cyberdyne/Project_Setup/* /scratch/ucgd/lustre/UCGD_Datahub/IRBs/IRB_00000000/Staging/
```
* Delete:
    * Ubox project
    * Mosaic project
    * Fabric project
    * Run the following commands to remove test project from the ucgddb:

```
$> psql -h remington.genetics.utah.edu -U UCGD-pepipeline -d ucgd -c "delete from projects where project = 'T101-101-NEO-Cyberdyne';"
$> psql -h remington.genetics.utah.edu -U UCGD-pepipeline -d ucgd -c "delete from samples where project = 'T101-101-NEO-Cyberdyne';"
```

## Notes

Project specific information can be found by running the following helper script:

`UCGDInspect projects -p [PROJECT]`
