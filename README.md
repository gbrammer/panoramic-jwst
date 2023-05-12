![Alt text](./PanoramicLogo.png?raw=true "PanoramicLogo")


# panoramic-jwst

Notebooks for bookkeeping NIRCam data from the PANORAMIC JWST survey.

*Run notebooks locally or on a GitHub codespace!*

[GO 2514 visit status report](https://www.stsci.edu/cgi-bin/get-proposal-info?id=2514&observatory=JWST)

# Basic workflow

- `Notebooks/JWST-associations.ipynb`
  - Query MAST for new data for a particular program ID (PANORAMIC = GO 2514)
  - Parse into grizli "associations" split by filter/detector
  - Send new associations to the grizli database
  - Optionally set non-default parameters for assoc pipeline processing
- `Notebooks/make-astrometry-catalog.ipynb`
  - Quicklook processing of a LW filter
  - Generate an astrometric catalog and send to the ``astrometry_reference`` table in the database
- Notebook TBD
  - *After* astrometric catalogs have been updated, set ``status=0`` for new associations
  - Launch AWS EC2 insstances to run preprocessing (alignment, etc.) on the ``status=0`` visits
- `Notebooks/panoramic-mosaics.ipynb`
  - Make 20/40 mas SW/LW mosaics of each PANORAMIC field
  - Send results to AWS/S3
  
## Mosaic products

https://s3.amazonaws.com/grizli-panoramic/mosaics/mosaics.html
