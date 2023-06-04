LCI_Analysis
================
Finn Lobnow
, Last edited on 31 May 2023

- <a href="#cell-lines" id="toc-cell-lines">CELL LINES</a>
- <a href="#parameters" id="toc-parameters">PARAMETERS</a>

## CELL LINES

- BDLD_27H-MyD88-TIR-TRAF6-BD-GFP TRAF6
- BDLD_50H-MyD88-TIR-TRAF6-BD-GFP TRAF6
- BDLD_14H-MyD88-TIR-TRAF6-BD-GFP TRAF6
- 3E10

## PARAMETERS

| parameter              | value | comments       |
|:-----------------------|------:|:---------------|
| tiff_compression_level |     5 | out of 10      |
| cell_diameter          |    25 | px, odd number |
| puncta_diameter        |     5 | px, odd number |

constants

| image                                    | exposure |
|:-----------------------------------------|:---------|
| 20220530 Darkfield 60 ms binned 001.tif  | 60 ms    |
| 20220530 Darkfield 100 ms binned 002.tif | 100 ms   |
| 20220530 Darkfield 200 ms binned 005.tif | 200 ms   |
| 20220530 Darkfield 300 ms binned 003.tif | 300 ms   |
| 20220530 Darkfield 400 ms binned 004.tif | 400 ms   |

dark_frames

| contains    | path                                                                                  |
|:------------|:--------------------------------------------------------------------------------------|
| input       | /Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1\_20230519/Input      |
| processing  | /Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1\_20230519/Processing |
| output      | /Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/batch_1\_20230519/Output     |
| dark_frames | /Volumes/TAYLOR-LAB/Finn/new_pipeline/pending_processing/dark_frames                  |
| ImageJ      | /Volumes/TAYLOR-LAB/Finn/Fiji.app/ImageJ-linux64                                      |

directories

| value       |
|:------------|
| Brightfield |
| IL-1        |

exclusion_channels

| image                                            | cohort                                | segment_with | ligand    | ligand_density | trackmate_max_link_distance | trackmate_threshold | trackmate_frame_gap | T Cy5 Finn protein_name | T GFP Finn protein_name | T RFP Cy3 Finn protein_name | WideField FIn protein_name |
|:-------------------------------------------------|:--------------------------------------|:-------------|:----------|---------------:|----------------------------:|--------------------:|--------------------:|:------------------------|:------------------------|:----------------------------|:---------------------------|
| 20230519 3nM_cl311-BDLD27H_TRAF6_MyD88 001.nd2   | BDLD_27H-MyD88-TIR-TRAF6-BD-GFP TRAF6 | MyD88        | 3 nM IL-1 |             20 |                         5.0 |                 2.0 |                   2 | IL-1                    | MyD88                   | TRAF6                       | Brightfield                |
| 20230519 3nM_cl313-BDLD50H_TRAF6_MyD88 001.nd2   | BDLD_50H-MyD88-TIR-TRAF6-BD-GFP TRAF6 | MyD88        | 3 nM IL-1 |             20 |                         5.0 |                 2.0 |                   2 | IL-1                    | MyD88                   | TRAF6                       | Brightfield                |
| 20230519 3nM_cl315-BDLD14H_TRAF6_MyD88 001.nd2   | BDLD_14H-MyD88-TIR-TRAF6-BD-GFP TRAF6 | MyD88        | 3 nM IL-1 |             20 |                         5.0 |                 2.0 |                   2 | IL-1                    | MyD88                   | TRAF6                       | Brightfield                |
| 20230519 GFP calibration_10pct_60ms 001.nd2      | Calibrations                          | GFP          |           |             NA |                         2.5 |                 2.0 |                   5 | IL-1                    | GFP                     | mScarlet                    | Brightfield                |
| 20230519 mScarlet calibration_2pct_300ms 001.nd2 | Calibrations                          | mScarlet     |           |             NA |                         2.5 |                 2.5 |                   5 | IL-1                    | GFP                     | mScarlet                    | Brightfield                |

images
