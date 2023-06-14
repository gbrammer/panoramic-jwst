"""
Scripts to make a FITSMap webpage and catalog overlays
"""

import os
import numpy as np
from grizli.aws import db
from grizli import utils

S3_MAP_PREFIX = 'grizli-v2/ClusterTiles/Map'

# Vizier catalogs
VIZ_CATALOGS = {'DESI-N (Duncan+22)':('VII/292/north','olive'),
                'DESI-S (Duncan+22)':('VII/292/south','olive'),
                'GAIA DR3': ('I/355/gaiadr3', 'lightblue'),
                'XMM-DR11 (Traulson+22)': ('IX/66/xmm411st', 'pink'),
                'XMM-DR12 (Webb+23)':('IIX/68/xmm4d12s', 'pink'),
                'SDSS DR16Q (Lyke+20)': ('VII/289/superset', 'magenta')
                #'Chandra (Evans+2019)': ('IX/57/csc2master', 'pink'), # Doesn't work?
               }


def run_make_fitsmap_html(field, viz_catalogs=VIZ_CATALOGS):
    """
    Run all steps to make the FITSMap HTML file and the catalog overlays
    """

    print(f'## {field}')

    html_file = f'{field}_map.html'

    f_tiles = db.SQL(f"""select * from combined_tiles_filters where field = '{field}'""")
    un = utils.Unique(f_tiles['filter'], verbose=True)

    ctiles = db.SQL(f"""select * from combined_tiles
    where field = '{field}'
    and tile = '09.09'
    """)

    crval1 = ctiles['crval1'][0]
    crval2 = ctiles['crval2'][0]

    make_html_file(field, crval1, crval2, f_tiles)

    make_tile_overlay(field)

    if viz_catalogs is not None:
        make_vizier_overlay(field, viz_catalogs=viz_catalogs, radius=10, with_desi=True)


def make_html_file(field, crval1, crval2, f_tiles):
    """
    Make FITSMap HTML file for a given field
    """

    html_file = f'{field}_map.html'

    HEAD = """<!DOCTYPE html>
<html>
<head>
   <title>FitsMap</title>
   <meta charset="utf-8" />
   <meta name="viewport" content="width=device-width, initial-scale=1.0">
   <link rel="shortcut icon" type="image/x-icon" href="docs/images/favicon.ico" />
   <link rel="stylesheet" href="https://unpkg.com/leaflet@1.3.4/dist/leaflet.css" integrity="sha512-puBpdR0798OZvTTbP4A8Ix/l+A4dHDD0DGqYW6RQ+9jxkRFclaxxQb/SJAWZfWAkuyeQUytO7+7N4QKrDh+drA==" crossorigin=""/>

   <script src='https://unpkg.com/leaflet@1.3.4/dist/leaflet.js' integrity='sha512-nMMmRyTVoLYqjP9hrbed9S+FzjZHW5gY1TWCHA5ckwXZBadntCNs8kEqAWdrb9O7rxbCaA4lKTIWjDXZxflOcA==' crossorigin=''></script>
   <script type='text/javascript' src='https://code.jquery.com/jquery-3.6.0.slim.min.js' charset='utf-8'></script>
    """

    HEAD += f"""
<script src='{field}_tiles.js'> </script>
<script src='{field}_vizier.js'> </script>
"""

    HEAD += """
   <style>
       html, body {
       height: 100%;
       margin: 0;
       }
       #map {
           width: 100%;
           height: 100%;
       }
       
       .leaflet-popup-content-wrapper {
           width: 530px;
       }
       .leaflet-popup-content {
           width: 520px;
       }
       
       #query {
           position: fixed;
           bottom: 0px;
           left: 0px;
           padding: 2px 2px 2px 2px;
           margin: 0;
           font-size: 8pt;
           font-family: Verdana;
           color: #333;
           z-index: 1000;
           background: #fff;
           background: rgba(255, 255, 255, 0.7);
       }
       
       #query a:link {
           text-decoration: none;
           color: #0078A8;
       }
       
       #query a:hover {
           text-decoration: underline;
           color: #0078A8;
       }

       input {
           font-size: 8pt;
       }
       
       #coordbox {
           position: fixed;
           bottom: 16px;
           left: 0px;
           padding: 2px 2px 2px 2px;
           margin: 0;
           font-size: 8pt;
           font-family: Verdana;
           color: #333;
           z-index: 1000;
           background: #fff;
           background: rgba(255, 255, 255, 0.7);
       }
       
   </style>
</head>
<body>
   <div id="map"></div>
   <div id="coordbox"> <input type="text" placeholder="RA Dec" id="boxtext" class="input" onkeydown="panFromBox(map)"/> </div>
   <div id="query"> - | - | - </div>
   
   <script>
    L.CRS.FitsMap = L.extend({}, L.CRS.Simple, {
        transformation: new L.Transformation(1/16, 0, -1/16, 256)
    });
    
    // var overlays = {}; // initialized in tiles.js
    """

    FOOTER = f"""
    var layerControl = L.control.layers(baseLayers, overlays);

    // WCS parameters
    // {field}
    var crpix = [1152.0, 1152.0]; 
    var crval = [{crval1:.6f}, {crval2:.6f}]; 
    var cdmatrix = [[-2.222222e-05, 0.000000e+00],
                    [0.000000e+00, 2.222222e-05]];
    """

    FOOTER += """
    var map = L.map("map", {
        crs: L.CRS.FitsMap,
        zoom: 2,
        minZoom: 0,
        preferCanvas: true,
    """

    FOOTER += f"""layers: {field.replace('-','_').lower()}_080_ncrgb,"""

    FOOTER += """
    });
    layerControl.addTo(map);
    
    var markers = null;

    urlParam = function(name){
        // Parse parameters from window.location, 
        // e.g., .../index.html?zoom=8
        // urlParam(zoom) = 8
        var results = new RegExp('[\?&]' + name + '=([^&#]*)').exec(window.location.href);
        if (results==null){
           return null;
        }
        else{
           return decodeURI(results[1]) || 0;
        }
    }

    pixToSky = function(xy){
        // Convert from zero-index pixel to sky coordinate assuming 
        // simple North-up WCS
        if (xy.hasOwnProperty('lng')){
            var dx = xy.lng - crpix[0] + 1;
            var dy = xy.lat - crpix[1] + 1;        
        } else {
            var dx = xy[0] - crpix[0] + 1;
            var dy = xy[1] - crpix[1] + 1;
        }
        var dra = dx * cdmatrix[0][0];
        var ddec = dy * cdmatrix[1][1];
        var ra = crval[0] + dra / Math.cos(crval[1]/180*3.14159);
        var dec = crval[1] + ddec;
        return [ra, dec];
    }

    skyToPix = function(rd){
        // Convert from sky to zero-index pixel coordinate assuming 
        // simple North-up WCS
        var dx = (rd[0] - crval[0]) * Math.cos(crval[1]/180*3.14159);
        var dy = (rd[1] - crval[1]);
        var x = crpix[0] - 1 + dx / cdmatrix[0][0];
        var y = crpix[1] - 1 + dy / cdmatrix[1][1];
        return [x,y];
    }

    skyToLatLng = function(rd){
        // Convert from sky to Leaflet.latLng coordinate assuming 
        // simple North-up WCS    
        var xy = skyToPix(rd);
        return L.latLng(xy[1], xy[0]);
    }

    setZoomFromUrl = function(map){
        // Read "zoom" from url bar
        var zoom = urlParam('zoom');
        if (zoom !== null){
            // console.log("zoom to " + zoom);
            map.setZoom(zoom);
        } else {
            map.setZoom(3);
        }
    }
    
    panToSky = function(rd, map){
        // Pan map to celestial coordinates
        var ll = skyToLatLng(rd)
        map.panTo(ll, map.getZoom());
        // console.log("pan to: " + rd + " / ll: " + ll.lng + ',' + ll.lat);
    }

    parseCoord = function(rd, hours) {
        // Parse sexagesimal coordinates if ':' found in rd
        if (rd.includes(':')){
            var dms = rd.split(':')
            
            var deg = dms[0]*1; 
            if (deg < 0) {
                var sign = -1;
            } else {
                var sign = 1;
            }
            deg += sign*dms[1]/60. + sign*dms[2]/3600.;
            if (hours > 0) {
                deg *= 360/24.
            }
        } else {
            var deg = rd;
        }
        return deg
    }

    panFromUrl = function(map){
        // Pan map based on ra/dec/[zoom] variables in location bar
        var ra = urlParam('ra');
        var dec = urlParam('dec');
        var coord = urlParam('coord');
        var coords = urlParam('coords');
        
        if ((ra !== null) & (dec !== null)) {
            panToSky([parseCoord(ra, 1), parseCoord(dec, 0)], map);
        } else if ((coord !== null)) {
            var rd = coord.split(',');
            if ((rd.length == 1)) {
                rd = coord.split(' ');
            }
            // console.log(rd);
            panToSky([parseCoord(rd[0], 1), parseCoord(rd[1], 0)], map)
        } else if ((coords !== null)) {
            var rd = coords.split(',');
            if ((rd.length == 1)) {
                rd = coords.split(' ');
            }
            // console.log(rd);
            panToSky([parseCoord(rd[0], 1), parseCoord(rd[1], 0)], map)
        } else {
            // Pan to crval
            panToSky(crval, map);
        }
        setZoomFromUrl(map);
    }
    
    panFromBox = function(map){
        // Pan map based on ra/dec/[zoom] variables in location bar
        var coord = $('#boxtext').val();
        var rd = coord.split(',');
        if ((rd.length == 1)) {
            rd = coord.split(' ');
        }
        panToSky([parseCoord(rd[0], 1), parseCoord(rd[1], 0)], map);
    }
    
    updateLocationBar = function(){
        var rd = pixToSky(map.getCenter());
        var params = 'coord=' + rd[0].toFixed(7);
        params += ',' + rd[1].toFixed(7);
        params += '&zoom=' + map.getZoom();
        var param_url = window.location.href.split('?')[0] + '?' + params;
        window.history.pushState('', '', param_url);

        // Query bar
        var legacyurl = 'https://www.legacysurvey.org/viewer?layer=ls-dr9&zoom=17';
        legacyurl += '&ra=' + rd[0].toFixed(7);
        legacyurl += '&dec=' + rd[1].toFixed(7);
        
        var query_html = '<a href="'+legacyurl+'">LegacySurvey</a>';
        
        var cdsurl = 'http://vizier.u-strasbg.fr/viz-bin/VizieR?&-c.rs=1&-c=';
        cdsurl += rd[0].toFixed(7);
        cdsurl += ',' + rd[1].toFixed(7);
        query_html += ' | <a href="'+cdsurl+'">CDS</a>';
        
        var esourl = 'https://archive.eso.org/scienceportal/home?pos=';
        esourl += rd[0].toFixed(7);
        esourl += ',' + rd[1].toFixed(7);
        esourl += '&r=0.02&dp_type=IMAGE,CUBE';
        query_html += ' | <a href="'+esourl+'">ESO</a>';
        
        var mfspectra = 'https://grizli-cutout.herokuapp.com/mosfire?mode=table&sep=5';
        mfspectra += '&ra=' + rd[0].toFixed(7);
        mfspectra += '&dec=' + rd[1].toFixed(7);
        query_html += ' | <a href="'+mfspectra+'">MOSFIRE</a>';
        
        var cutout = 'https://grizli-cutout.herokuapp.com/thumb?all_filters=True&size=4&scl=1.0&asinh=True&filters=f115w-clear,f277w-clear,f444w-clear&rgb_scl=1.5,0.74,1.3&pl=2';
        cutout += '&ra=' + rd[0].toFixed(7);
        cutout += '&dec=' + rd[1].toFixed(7);
        query_html += ' | <a href="'+cutout+'">NIRCam</a>';
    """

    ### Layers by filter
    un = utils.Unique(f_tiles['filter'], verbose=False)

    all_filters = ','.join([v.lower() for v in un.values])
    # if 'f814w' not in all_filters:
    #     all_filters = 'f814w,' + all_filters

    FOOTER += f"""
        var cut2 = 'https://grizli-cutout.herokuapp.com/thumb?all_filters=True&size=4&scl=1.0&asinh=False&filters={all_filters}&rgb_scl=1.1,1.05,1.0&pl=2';
        cut2 += '&ra=' + rd[0].toFixed(7);
        cut2 += '&dec=' + rd[1].toFixed(7);
        query_html += ' | <a href="'+cut2+'">All filters</a>';
    """

    FOOTER += """
        // var miricut = 'https://grizli-cutout.herokuapp.com/thumb?all_filters=True&size=4&scl=1&asinh=True&filters=f150w-clear,f444w-clear,f770w&rgb_scl=0.8,1.2,1.0&pl=2.0';
        // miricut += '&ra=' + rd[0].toFixed(7);
        // miricut += '&dec=' + rd[1].toFixed(7);
        // query_html += ' | <a href="'+miricut+'">NRC+MIR</a>';

        $('#query').html(query_html);
    }
    
    clearOverlays = function(){
        if (markers !== null){
            for(i = 0; i < markers.length; i++) {
                overlays[labels[i]].remove();
                console.log('clear overlays: ' + labels[i])
            }
        }
    }

    $( document ).ready(function() {    
        // Set pan from URL bar, if ra/dec/zoom specified
        panFromUrl(map);
        
        // Clear any catalog overlays
        clearOverlays();
        
        // Event handlers to update location bar when view changes
        map.on('moveend', updateLocationBar);
        map.on('zoomend', updateLocationBar);

    });    
    
   </script>
</body>
</html>
    """

    f_low = field.replace('-','_').lower()

    with open(html_file,'w') as fp:
        fp.write(HEAD)
        
        field_filters = un.values + ['swrgb','lwrgb','ncrgb']
        
        for f in field_filters:
            fi = f.lower()
            line = f"""var {f_low}_080_{fi.split('-clear')[0]} =  L.tileLayer("{f_low}_080_{fi}"""
            line += """/{z}/{y}/{x}.png",{ attribution:"""
            line += f'"{f}"'
            line += """, minZoom: 0, maxZoom: 9, minNativeZoom: 1, maxNativeZoom: 4,});"""
            fp.write(line + '\n')
        
        fp.write("""
        var baseLayers = {\n""")
        
        for f in field_filters:
            fi = f.lower()
            fp.write(f"""        "{f}": {f_low}_080_{fi.split('-clear')[0]},\n""")
        
        fp.write("    };\n")
        
        fp.write(FOOTER)
        
    os.system(f'aws s3 cp {html_file} s3://{S3_MAP_PREFIX}/{field}/index.html --acl public-read')
    print(f"https://s3.amazonaws.com/{S3_MAP_PREFIX}/{field}/index.html")


def get_tile_wcs(field, ref='09.09'):
    """
    Get WCS of reference tile
    """
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs

    ref_tile = db.SQL(f"""select * from combined_tiles where field = '{field}' AND tile = '{ref}'""")

    h = pyfits.Header()
    for k in ref_tile.colnames:
        h[k[:8]] = ref_tile[k][0]
        
    wcs = pywcs.WCS(h)

    return wcs


def make_tile_overlay(field):
    """
    Make a javascript overlay with links to the tile FITS files.
    
    """

    tiles = db.SQL(f"""select c.field, c.tile, c.filter, t.footprint from combined_tiles_filters c, combined_tiles t
    where (c.field = t.field) AND (c.tile = t.tile) AND c.field = '{field}'
    """)
        
    wcs = get_tile_wcs(field, ref='09.09')

    key = ['{field}-080-{tile}'.format(**row) for row in tiles]
    un = utils.Unique(key, verbose=False)
    ix = un.unique_index()

    lix = un.list_indices

    with open(f'{field}_tiles.js','w') as fp:
        fp.write('var overlays = {};\nvar tiles = [];\n')
        
        for i, (v, f) in enumerate(zip(un.values, tiles['footprint'][ix])):
            filters = tiles['filter'][lix[i]].tolist()
            
            so = np.argsort(filters)
            filters = [filters[i] for i in so]
            
            print(v, ' '.join(filters))
            
            sr = utils.SRegion(f)
            xy = wcs.all_world2pix(*sr.xy[0].T, 0)
            xyl = np.array(xy[::-1]).T.tolist()
            xyls = ','.join(['[{0:.1f},{1:.1f}]'.format(*r) for r in xyl])
            
            tt = f"""bindTooltip('{v}', {{ direction: 'auto'}})"""
            path = f"https://s3.amazonaws.com/{S3_MAP_PREFIX.replace('/Map','')}/{field}"
            pop = f"bindPopup('<h3> Tile {v.split('-')[-1]} </h3> "
            for f in filters:
                if f in ['F110W','F125W','F140W','F160W','F770W','F1800W']:
                    file = f'{v}-{f.lower()}_drz_sci.fits.gz'
                else:
                    file = f'{v}-{f.lower()}_drc_sci.fits.gz'
                
                pop += f'<a href="{path}/{file}" /> {file} </a> '
                pop += f"/ <a href=\"{path}/{file.replace('sci','wht')}\" /> wht </a> <br>"

            color = 'white'
            
            if 'F444W-CLEAR' in filters:
                color = 'cornsilk'
            if 'F770W' in filters:
                color = 'crimson'
                
            pop += "')"
            fp.write(f"""tiles.push(L.polygon([ {xyls} ], {{color: '{color}', weight:2, opacity:0.2, fill:true}}).{tt}.{pop});\n""")

            # if i > 5:
            #     break
        
        fp.write("overlays['FITS tiles'] = L.layerGroup(tiles);")
        
    os.system(f'aws s3 cp {field}_tiles.js s3://{S3_MAP_PREFIX}/{field}/ --acl public-read')


def make_vizier_overlay(field, viz_catalogs=VIZ_CATALOGS, ref_tile='09.09', radius=10, with_desi=True):
    """
    Make a JavaScript overlay with some Vizier catalog queries
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex

    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs

    from grizli import utils
    import grizli.catalog

    v = 'J/A+A/602/A1/table1'

    output_file = f'{field}_vizier.js'
    if os.path.exists(output_file):
        os.remove(output_file)

    has_catalog = []

    wcs = get_tile_wcs(field, ref=ref_tile)

    r0, d0 = wcs.calc_footprint().mean(axis=0)
    ra_ref, dec_ref = r0, d0

    print(f'# {field}')

    for vname in viz_catalogs:
        if (field == 'gds') & ('DESI' in vname):
            continue
            
        if (vname == 'DESI-S (Duncan+22)') & ('DESI-N (Duncan+22)' in has_catalog):
            print('    - Skip ', vname)
            continue

        v, vcolor = viz_catalogs[vname]
        
        if 'GAIA' in vname:
            try:
                # GAIA evaluated at MJD=60000 (2023)
                vcat = grizli.catalog.get_gaia_DR2_vizier(ra=ra_ref, dec=dec_ref,
                                                    radius=radius, mjd=60000, clean_mjd=False)
            except:
                print(f"Query {vname}: {v} failed")
                continue
        else:
            try:
                vcat = grizli.catalog.query_tap_catalog(ra=r0, dec=d0,
                                        radius=radius, vizier=False,
                                        tap_url="http://tapvizier.u-strasbg.fr/TAPVizieR/tap/",
                                        #tap_url="http://TAPVizieR.u-strasbg.fr/TAPVizieR/tap/",
                                        rd_colnames = ['RAJ2000', 'DEJ2000'],
                                        db=f'"{v}"', verbose=False)
            except:
                print(f"    * Query {vname}: {v} failed")
                continue
            
        print('    ', vname, v, len(vcat))

        if len(vcat) == 0:
            continue
        
        has_catalog.append(vname)
            
        cm = plt.cm.summer
        
        popups = None

        if 'zphot' in vcat.colnames:
            vcat['comment'] = ['id = {id} <br> z<sub>phot</sub> = {zphot:.2f}'.format(vname=vname, **row) for row in vcat]
        elif 'DR3Name' in vcat.colnames:
            vcat['comment'] = ['{DR3Name} <br> pm = {pmRA:.2f}, {pmDE:.2f} mas'.format(**row) for row in vcat]
        elif 'pmra' in vcat.colnames:
            vcat['comment'] = ['pm = {pmra:.2f}, {pmdec:.2f} mas'.format(**row) for row in vcat]
            vcat['ra_time'][~np.isfinite(vcat['dec_time'])] = vcat['ra'][~np.isfinite(vcat['dec_time'])]
            vcat['dec_time'][~np.isfinite(vcat['dec_time'])] = vcat['dec'][~np.isfinite(vcat['dec_time'])]
            vcat['ra'] = vcat['ra_time']
            vcat['dec'] = vcat['dec_time']
        
        elif 'Lyke' in vname:
            # quasars
            vcat['comment'] = ['SDSS = {SDSS} &nbsp;&nbsp; z = {z:.3f}'.format(vname=vname, **row) for row in vcat]
            popups = ['SDSS = {SDSS} &nbsp;&nbsp; z = {z:.3f} &nbsp;&nbsp; <a href="http://dr16.sdss.org/optical/spectrum/view?plateid={Plate}&mjd={MJD}&fiberid={Fiber}" /> Spectrum </a>'.format(vname=vname, **row)
                      for row in vcat]

        elif '2CXO' in vcat.colnames:
            vcat['comment'] = ['Chandra Source Catalog <br> 2CXO = {2CXO}'.format(**row) for row in vcat]
            
        elif 'MUSE-Wide' in vcat.colnames:
            vcat['comment'] = ['{MUSE-Wide} {Field}  z={z:.4f} ({OLines}, {LeadLine} SN={SN:.1f})'.format(**row)
                            for row in vcat]
        elif 'MUSE' in vcat.colnames:
            vcat['comment'] = ['{MUSE} {r_zMuse}  zMUSE={zMuse:.4f}'.format(**row)
                            for row in vcat]

        elif 'DLyAFit' in vcat.colnames:
            vcat['comment'] = ['{ID} {DataSet}  zMUSE={z:.4f}'.format(**row)
                            for row in vcat]

        else:
            vcat['comment'] = vname
        
        sizes = None
        if ('XMM' in vname):
            for c in ['ePos','srcML']:
                if c in vcat.colnames:
                    sizes = vcat[c]
                    print(f'Use {c} for source sizes')
                    break

        #col = cm(np.clip((vcat['SNR']-5)/8, 0,1)*0.7+0.2)
        
        key = vname.split()[0].replace('-','_').lower()
        
        with open(output_file,'a') as fp:
            fp.write(f'var {key} = [];\n')

            olay = []
            si = 0.5/0.1 # R=0.5"
            col = vcolor
            
            for i, (ri, di, comment) in enumerate(zip(vcat['ra'], vcat['dec'], vcat['comment'])):

                if sizes is None:
                    si = 0.5/0.1
                else:
                    si = sizes[i]/0.1

                sr = utils.SRegion(f'CIRCLE({ri},{di},{si}")', wrap=False, ncircle=16)
                xy = wcs.all_world2pix(*sr.xy[0].T, 0)
                xyl = np.array(xy[::-1]).T.tolist()
                xyls = ','.join(['[{0:.1f},{1:.1f}]'.format(*r) for r in xyl])

                xpi, ypi = wcs.all_world2pix([ri], [di], 0)

                # olay.append(f'[{xyls}]')

                if sizes is None:
                    marker = f"""L.circleMarker([{ypi[0]:.1f},{xpi[0]:.1f}], {{radius:8,color:'{col}',weight:2,opacity:0.8,fill:false}}).bindTooltip('{comment}', {{direction:'auto'}})"""
                else:
                    marker = f"""L.polygon([ {xyls} ], {{color: '{col}', weight:2, opacity:0.8, fill:false}}).bindTooltip('{comment}', {{ direction: 'auto'}})"""
                
                if popups is not None:
                    marker += f""".bindPopup('{popups[i]}', {{ direction: 'auto'}})"""

                fp.write(f"""{key}.push({marker});\n""")

            fp.write(f"""overlays['{vname}'] = L.layerGroup({key});""")
        
        os.system(f'aws s3 cp {output_file} s3://{S3_MAP_PREFIX}/{field}/ --acl public-read --quiet')

    if with_desi:
        _ = make_desi_edr_layer(field, ref_tile=ref_tile)

    print(f"     https://s3.amazonaws.com/{S3_MAP_PREFIX}/{field}/index.html")


def make_desi_edr_layer(field, ref_tile='09.09'):
    """
    """
    import grizli.catalog

    output_file = f'{field}_vizier.js'

    wcs = get_tile_wcs(field, ref=ref_tile)
    r0, d0 = wcs.calc_footprint().mean(axis=0)
    ra_ref, dec_ref = r0, d0

    # DESI
    edr = grizli.catalog.query_tap_catalog(ra=ra_ref, dec=dec_ref, radius=30,
                                       tap_url='https://datalab.noirlab.edu/tap',
                                       db='desi_edr.zpix',
                                       rd_colnames=['mean_fiber_ra', 'mean_fiber_dec'],
                                       verbose=False,
                                )
    
    if len(edr) == 0:
        return False
    
    else:
        print(f'     DESI EDR: {len(edr)} sources')

    key = 'desi_edr'
    
    edr['comment'] = ['targetid={targetid} <br> {spectype} z={z:.4f} &plusmn; {zerr:.4f}'.format(**row) for row in edr]
    edr['popup'] = ['{comment} <br> <a href="https://www.legacysurvey.org/viewer/desi-spectrum/edr/targetid{targetid}"/> Spectrum </a>'.format(**row)
                     for row in edr]
    
    sizes = None

    with open(output_file,'a') as fp:
        fp.write(f'var {key} = [];\n')

        olay = []
        si = 0.5/0.1 # R=0.5"
        
        for i in range(len(edr)):
            # , (ri, di, comment) in enumerate(zip(vcat['ra'], vcat['dec'], vcat['comment'])):
            ri = edr['ra'][i]
            di = edr['dec'][i]
            comment = edr['comment'][i]
            popup = edr['popup'][i]

            if sizes is None:
                si = 0.5/0.1
            else:
                si = sizes[i]/0.1

            sr = utils.SRegion(f'CIRCLE({ri},{di},{si}")', wrap=False, ncircle=16)
            xy = wcs.all_world2pix(*sr.xy[0].T, 0)
            xyl = np.array(xy[::-1]).T.tolist()
            xyls = ','.join([f'[{r[0]:.1f},{r[1]:.1f}]' for r in xyl])

            xpi, ypi = wcs.all_world2pix([ri], [di], 0)

            # olay.append(f'[{xyls}]')

            if edr['spectype'][i] == 'STAR':
                col = 'lightblue'
            elif edr['spectype'][i] == 'GALAXY':
                col = 'yellow'
            elif edr['spectype'][i].upper() in ['QSO','QUASAR']:
                col = 'magenta'
            else:
                col = 'white'

            if sizes is None:
                marker = f"""L.circleMarker([{ypi[0]:.1f},{xpi[0]:.1f}], {{radius:8,color:'{col}',weight:2,opacity:0.8,fill:false}}).bindTooltip('{comment}', {{direction:'auto'}})"""
            else:
                marker = f"""L.polygon([ {xyls} ], {{color: '{col}', weight:2, opacity:0.8, fill:false}}).bindTooltip('{comment}', {{ direction: 'auto'}})"""
            
            marker += f""".bindPopup('{popup}', {{ direction: 'auto'}})"""

            fp.write(f"""{key}.push({marker});\n""")

        fp.write(f"""overlays['DESI EDR'] = L.layerGroup({key});""")
    
    os.system(f'aws s3 cp {output_file} s3://{S3_MAP_PREFIX}/{field}/ --acl public-read --quiet')
    

