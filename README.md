![Alt text](./PanoramicLogo.png?raw=true "PanoramicLogo")


# panoramic-jwst

Notebooks for bookkeeping NIRCam data from the PANORAMIC JWST survey.

*Run notebooks locally or on a GitHub codespace!*

[GO 2514 visit status report](https://www.stsci.edu/cgi-bin/get-proposal-info?id=2514&observatory=JWST)

*NB: If you're not part of the collaboration, the steps below that require write access to the grizli database won't work.  However, the general processing scripts should work on a standalone machine with a standard grizli installation, e.g., `pip install grizli[aws,jwst]`.*

# Basic workflow

1. `Notebooks/JWST-associations.ipynb`
  - Query MAST for new data for a particular program ID (PANORAMIC = GO 2514)
  - Parse into grizli "associations" split by filter/detector
  - Send new associations to the grizli database
  - **In detail**:
    1. `Run All` in `Notebooks/JWST-associations.ipynb` will generate the associations and stop at a breakpoint labeled *Break here....*
    2. Check that the identified associations are complete, i.e., all filter+detector combinations have the same number of files.
    3. Check at *Check for rows that are already in the database...* if there were new associations found.  If so, the message will print "New: **N**" with **N** > 0.  
    4. To send the new associations to the database, flip the "False" to "True" in the first `if/else` clause at *# Send associations to the database table*.
    5. *Flip that clause back to False after running to avoid sending duplicate association definitions.*. It won't be a problem if you "Run All" the notebook again since they'll be identified as already having been sent, but if you run the "send to db" line multiple times consecutively it could create duplicate rows in the association table.
2. `Notebooks/make-astrometry-catalog.ipynb`
  - Quicklook processing of a LW filter
  - Generate an astrometric catalog and send to the ``astrometry_reference`` table in the database.
  - The first time an association is run that doesn't overlap with any astrometric reference, the pipeline will try to pull astrometric reference files from some large ground-based database, e.g., PanSTARRS or NOAO LegacySurveys.  These alone seem to be reasonably well aligned to GAIA, but even their source density is fairly low within the JWST instrumental FOV.  So the procedure is to first align some LW filter to the external reference (LW because of high sensitivity and larger FOV).  Then make a catalog from that and send the sources to the `astrometry_reference` table in the databse, which will then be queried automatic when the other filters that overlap with that FoV are run.
  - **In detail**:
    1. `Run All` in `Notebooks/make-astrometry-catalog.ipynb`.  This will find PANORAMIC associations that haven't yet been processed and will run the F444W associations for them and generate the catalogs.
    2. Check the lines that print the contents of the `shifts.log` and `wcs.log` files.  If the associations were aligned 
      - The `shifts.log` files should indicate that at least a few dozen sources were identified in the different exposures and were aligned with RMS <~ 0.1 pixel
      - The `wcs.log` files will indicate the reference catalog to which the association was aligned.  E.g., `# radec: indef-02514-209-124.0-nrcb5-f444w-clear_ls_dr9.radec` was aligned to `ls_dr9` = LegacySurveys DR9.  Check that the number of alignment sources *N* is at least a few and that the RMS is of order < 0.5-1 pixel.
    3. If the alignment logs look OK and if there were a few hundred sources found in each derived alignment catalog (e.g., "Add 578 objects to astrometry_reference for j131432p2432_indef-f115w..."), flip `SEND_TO_DB = False` under **Process reference catalogs** to send the reference source lists to the `astrometry_reference` table in the database. 
    4. Finally, run the cell with `with open('astrometry_log.txt','a') ...` and commit the changes to `astrometry_log.txt` to the repository. 
3. Still in `Notebooks/make-astrometry-catalog.ipynb`
  - *After* astrometric catalogs have been updated, set ``status=0`` for new associations
  - Optionally set non-default parameters for assoc pipeline processing
    - Good for 1/f processing in PANORAMIC and other relatively shallow programs
  - Launch AWS EC2 insstances to run preprocessing (alignment, etc.) on the ``status=0`` visits
  - **In detail**:
    1. Run the few cells at `## Reset status and launch EC2` to set the `status=0` of the unprocessed associations.
    2. Run the cell to make the pipeline parameter YAML files, which will be uploaded to the necessary location on S3
    3. Flip the switch `If 0:` > `If 1:` and run at `## Launch EC2 instances`.
    4. This will launch _N_ EC2 instances that will each run through the processing of the `status=0` associations.  It takes some time for the instances to spin up, but you can then monitor their process by rerunning the cells at `# Check status of associations in the DB`.  At first there will be _N_ > 0 entries with `status=0` in the query there, and eventually you'll see them start to change to `status=1` as they're processed.  After a few minutes, all of the `status=0` associations should be done and the number you saw originally should now be under `status=2`, meaning "Done". There are two copies of the status query cell to make it easier to run them one at a time and see how the numbers change.
4. `Notebooks/panoramic-mosaics.ipynb`
  - Make 20/40 mas SW/LW mosaics of each PANORAMIC field
  - Send results to AWS/S3
  
## Mosaic products

https://s3.amazonaws.com/grizli-panoramic/mosaics/mosaics.html
