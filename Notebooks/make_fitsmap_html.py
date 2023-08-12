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

    eso_query_layer(field)
    
    mosfire_slits_layer(field)
    
    nirspec_slits_layer(field)


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
<script src='{field}_eso.js'> </script>
<script src='{field}_mosfire.js'> </script>
<script src='{field}_nirspec.js'> </script>
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
    
    if field == 'abell2744':
        ref = '08.08'
    elif field == 'uds':
        ref = '11.10'
    elif field == 'cos':
        ref = '16.16'
    elif field.startswith('egs'):
        ref = '10.14'
    
    print(f'WCS reference tile {ref} for field {field}')
    
    ref_tile = db.SQL(f"""select * from combined_tiles where field = '{field}' AND tile = '{ref}'""")

    h = pyfits.Header()
    for k in ref_tile.colnames:
        h[k[:8]] = ref_tile[k][0]
        
    wcs = pywcs.WCS(h)

    return wcs


def make_tile_overlay(field, ref_tile='09.09'):
    """
    Make a javascript overlay with links to the tile FITS files.
    
    """

    tiles = db.SQL(f"""select c.field, c.tile, c.filter, t.footprint from combined_tiles_filters c, combined_tiles t
    where (c.field = t.field) AND (c.tile = t.tile) AND c.field = '{field}'
    """)
    
    if field == 'abel2744':
        ref_tile = '08.08'
    
    wcs = get_tile_wcs(field, ref=ref_tile)

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
                    print(f'       - Use {c} for source sizes')
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
    

def eso_query_layer(field, upload=True):
    
    import grizli.catalog
    import astropy.io.fits as pyfits
    
    output_file = f'{field}_eso.js'

    wcs = get_tile_wcs(field)
    r0, d0 = wcs.calc_footprint().mean(axis=0)
    ra_ref, dec_ref = r0, d0
    
    dd = 30. # arcmin
    dx = dd/60/np.cos(d0/180*np.pi)
    dy = dd/60
    
    # coord = "AND s_ra > 149. AND s_ra < 151 AND s_dec > 1.4 AND s_dec < 3.0" # COSMOS

    coord = f"AND s_ra > {r0-dx} AND s_ra < {r0+dx}"
    coord += f" AND s_dec > {d0-dy} AND s_dec < {d0+dy}" 

    query = f"""SELECT * FROM ivoa.ObsCore
    WHERE (obs_collection='ALMA' OR obs_collection='MUSE' OR obs_collection='MUSE-DEEP')
    {coord}
    """
    #query = "SELECT * FROM ivoa.ObsCore WHERE obs_collection='MUSE-DEEP' OR obs_collection='MUSE'"

    MAXREC=200000
    qstr = query.replace("'","%27").replace(' ','+').replace('\n','+').replace('=','%3D')
    qstr = qstr.replace('>','%3E').replace('<','%3C')
    TAP_URL = "http://archive.eso.org/tap_obs/sync?REQUEST=doQuery&FORMAT=fits&LANG=ADQL&MAXREC={0}&QUERY={1}"
    SEND = TAP_URL.format(MAXREC, qstr)
    print(SEND)

    # wcs = pywcs.WCS(_h)

    imq = pyfits.open(SEND)
    alma = utils.GTable()
    for c in imq[1].data.columns:
        if len(imq[1].data[c.name]) > 0:
            alma[c.name] = imq[1].data[c.name]

    alma['obs_collection'] = [o.strip() for o in alma['obs_collection']]

    muse = alma['obs_collection'] == 'MUSE'
    alma['filter'][muse] = '0'

    muse = alma['obs_collection'] == 'MUSE-DEEP'
    alma['filter'][muse] = '1'

    keep = np.array([s.startswith('POLY') for s in alma['s_region']])
    alma = alma[keep]
    
    print(f'Full ESO query: {len(alma)}')
    if len(alma) == 0:
        return True
        
    with open(output_file,'w') as fp:
    
        fp.write('    // ALMA \n')
    
        un = utils.Unique(alma['filter'], verbose=False)

        for k, sel in un:
            if sel.sum() > 100000:
                continue
            
            if k == '0':
                key = f'MUSE: {sel.sum()}'
                color = 'magenta'
            elif k == '1':
                key = f'MUSE-DEEP: {sel.sum()}'
                color = 'purple'
            else:
                key = f'ALMA Band {k}: {sel.sum()}'
                color='orange'
            
            fp.write(f"    var alma_overlays_{k} = [];\n")

            print(key, sel.sum())
            regs = ' '.join([p for p in alma['s_region'][sel]])
            sr = utils.SRegion(regs)
        
            ix = np.where(sel)[0]
            area = sr.sky_area()
            center = sr.centroid

            for i, rd in enumerate(sr.xy):
                xy = wcs.all_world2pix(rd, 0)
            
                poly = ','.join([f'[{c[1]:.1f}, {c[0]:.1f}]' for c in xy])

                if xy.shape[0] > 16:
                    sri = utils.SRegion('CIRCLE({0},{1},{2}")'.format(*center[i], np.sqrt(area[i].value/np.pi)*60.), wrap=False, ncircle=16)
                    xyi = wcs.all_world2pix(*sri.xy[0].T, 0)
                    xyl = np.array(xyi[::-1]).T.tolist()
                    xyls = ','.join(['[{0:.1f},{1:.1f}]'.format(*r) for r in xyl])
                
                    if np.allclose(area[i], sri.sky_area()[0], atol=0.1):
                        poly = xyls
                        
                title = alma['obs_title'][ix][i].strip().replace("'",'')
                pi = alma['obs_creator_name'][ix][i].split(',')[0].strip().replace("'",'')
                targ = alma['target_name'][ix][i].strip().replace("'",'')
                url = alma['access_url'][ix][i].strip()
            
                if 'ID?ADP' in url:
                    adp = url.split('ID?ADP')[1]
                    url = f'https://archive.eso.org/dataset/ADP{adp}'
                elif 'almascience' in url:
                    url += f'&observationsSourceName={targ}'
                
                prop = alma['proposal_id'][ix][i].strip()
            
                pop = f"'{title} <br> {prop} ({pi}) : <a href=\"{url}\" > {targ} </a>'"
                
                if k == '0':
                    tt = f".bindTooltip('MUSE {targ}', {{ direction: 'auto'}})"
                elif k == '1':
                    tt = f".bindTooltip('MUSE-DEEP {targ}', {{ direction: 'auto'}})"
                else:
                    tt = f".bindTooltip('ALMA Band {k} {targ}', {{ direction: 'auto'}})"

                #pop += f" <img src=\"https://s3.amazonaws.com/mosfire-pipeline/Spectra/{mf['file'][i].replace('sp.fits','sp2d.png')}\" width=500px />"
                #pop += f" <br> <img src=\"https://s3.amazonaws.com/mosfire-pipeline/Spectra/{mf['file'][i].replace('sp.fits','sp_log1d.png')}\" width=500px />'"
                row = f"\n    alma_overlays_{k}.push(L.polygon([{poly}],"
                row += f" {{color: '{color}', weight:1, opacity:0.8, fill:false}}).bindPopup({pop}, {{minWidth: 500, maxWidth:520}}){tt});"
                fp.write(row)
        
            row = f"    overlays['{key}'] = L.layerGroup(alma_overlays_{k});\n"
            fp.write(row)
    
    if upload:
        os.system(f'aws s3 cp {output_file} s3://{S3_MAP_PREFIX}/{field}/ --acl public-read')


def mosfire_slits_layer(field, upload=True):
    """
    Query database for MOSFIRE spectra
    """
    output_file = f'{field}_mosfire.js'

    wcs = get_tile_wcs(field)
    
    r0, d0 = wcs.calc_footprint().mean(axis=0)
    ra_ref, dec_ref = r0, d0
    
    dd = 30. # arcmin
    dx = dd/60/np.cos(d0/180*np.pi)
    dy = dd/60
    
    # coord = "AND s_ra > 149. AND s_ra < 151 AND s_dec > 1.4 AND s_dec < 3.0" # COSMOS

    mf = db.SQL(f"""select file, slit_width, slit_length, skypa3, targoff, ra_slit, dec_slit, filter, datemask, target_name
    from mosfire_extractions
    where ra_slit > {r0-dx} AND ra_slit < {r0+dx}
    AND dec_slit > {d0-dy} AND dec_slit < {d0+dy}""")

    print('MOSFIRE N: ', len(mf))
    
    if len(mf) == 0:
        return True
        
    with open(output_file,'w') as fp:

        row = '    var mosfire_slits = [];\n    var mosfire_targets = [];'
        fp.write(row)

        for i in range(len(mf)):
            sx = np.array([-0.5, -0.5, 0.5, 0.5, -0.5])*mf['slit_width'][i]/0.08
            sy = (np.array([-0.5, 0.5, 0.5, -0.5, -0.5])*mf['slit_length'][i])/0.08
            syo = (np.array([-0.5, 0.5, 0.5, -0.5, -0.5])*mf['slit_width'][i] + mf['targoff'][i])/0.08 # length of sp2d.png cutout
            theta = -mf['skypa3'][i]/180*np.pi
            _mat = np.array([[np.cos(theta), -np.sin(theta)],
                             [np.sin(theta), np.cos(theta)]])

            srot = np.array([sx, sy]).T.dot(_mat)
            srot += np.squeeze(wcs.all_world2pix([mf['ra_slit'][i]], 
                                                 [mf['dec_slit'][i]], 1))
            trot = np.array([sx, syo]).T.dot(_mat)
            trot += np.squeeze(wcs.all_world2pix([mf['ra_slit'][i]], 
                                                 [mf['dec_slit'][i]], 1))
        
            # plt.plot(*srot.T)
        
            poly = ','.join([f'[{c[1]:6.1f}, {c[0]:6.1f}]' for c in srot])
            tpoly = ','.join([f'[{c[1]:6.1f}, {c[0]:6.1f}]' for c in trot])
        
            pop = f"'{mf['datemask'][i]} <b> {mf['target_name'][i]} </b> ({mf['filter'][i]}) <br> "
            pop += f" <img src=\"https://s3.amazonaws.com/mosfire-pipeline/Spectra/{mf['file'][i].replace('sp.fits','sp2d.png')}\" width=500px />"
            pop += f" <br> <img src=\"https://s3.amazonaws.com/mosfire-pipeline/Spectra/{mf['file'][i].replace('sp.fits','sp_log1d.png')}\" width=500px />'"
            row = f"\n    mosfire_slits.push(L.polygon([{poly}],"
            row += f" {{color: 'pink', weight:1, opacity:0.8, fill:false}}).bindPopup({pop}));"
            fp.write(row+'\n')
            row = f"\n    mosfire_targets.push(L.polygon([{tpoly}],"
            row += f" {{color: 'magenta', weight:1, opacity:0.8, fill:false}}).bindPopup({pop}));"
            #print(row)
            fp.write(row+'\n')
    
        if len(mf) > 0:
            row = "    overlays['MOSFIRE slits'] = L.layerGroup(mosfire_slits);\n"
            fp.write(row)
            row = "    overlays['MOSFIRE targets'] = L.layerGroup(mosfire_targets);\n"
            fp.write(row)
    
    if upload:
        os.system(f'aws s3 cp {output_file} s3://{S3_MAP_PREFIX}/{field}/ --acl public-read')


def nirspec_slits_layer(field, upload=True):
    """
    Query database for NIRSpec slits
    """
    output_file = f'{field}_nirspec.js'

    wcs = get_tile_wcs(field)
    
    r0, d0 = wcs.calc_footprint().mean(axis=0)
    ra_ref, dec_ref = r0, d0
    
    dd = 30. # arcmin
    dx = dd/60/np.cos(d0/180*np.pi)
    dy = dd/60
    
    slits = db.SQL(f"""select program, msametfl, msametid, slitlet_id, grating, filter,
    ra, dec, source_id, is_source, footprint
    from nirspec_slits
    where patt_num = 1 
    AND polygon(circle(point({r0},{d0}),0.6)) @> point(ra, dec)
    """)
    
    if len(slits) == 0:
        return True
        
    keys = ['{program} {grating} {filter}'.format(**row).replace(' CLEAR','')
            for row in slits]
    
    un = utils.Unique(keys)
    
    rows = []
    
    rows.append("var sprop0 = {color: 'white', weight:1, opacity:0.8, fill:false};")
    rows.append("var sprop1 = {color: 'magenta', weight:1, opacity:0.8, fill:false};")
    
    for j, k in enumerate(un.values):
        #kl = 'slit_'+k.lower().replace(' ','_')
        kl = f'nrs_{j}'
        rows.append(f'var {kl} = []; // {k}')
        for row in slits[un[k]]:
            sr = utils.SRegion(row['footprint'])
            
            xy = wcs.all_world2pix(sr.xy[0], 0)
            srx = utils.SRegion(np.array(xy)[:,::-1], wrap=False)
            
            prop = 'sprop1' if row['is_source'] else 'sprop0'
            poly = srx.polystr(precision=1)[0].replace('(','[').replace(')',']')
            
            marker = f'L.polygon([{poly}],{prop})'
            
            if row['is_source']:
                marker += ".bindTooltip('{program} {grating} {filter} #{source_id}')".format(**row).replace(' CLEAR','')
                
            row = f"{kl}.push({marker});"
            #row += f" {{color: '{color}', weight:1, opacity:0.8, fill:false}}));"
            
            rows.append(row)
            
        rows.append(f"overlays['Slits {k}'] = L.layerGroup({kl});")
    
    
    # Extractions
    nre = db.SQL(f"""select root, file, ra, dec, grating, filter, SUBSTR(dataset,4,4) as program
    from nirspec_extractions
    WHERE polygon(circle(point({r0},{d0}),0.6)) @> point(ra, dec) and root not like 'uncover-5m%%'
    
    """)
    
    if len(nre) > 0:
        nre['grade'] = -1
        nre['z'] = -1.
    
        nrz = db.SQL(f"""select nz.root, nz.file, nz.z as zauto, nm.z as z, grade
        from nirspec_redshifts nz, nirspec_redshifts_manual nm, nirspec_extractions ne
        WHERE nz.file = nm.file AND ne.file = nz.file
        AND polygon(circle(point({r0},{d0}),0.6)) @> point(ra, dec) and nz.root not like 'uncover-5m%%'
    
        """)
    
        for i, row in enumerate(nrz):
            for k in ['grade','z']:
                nre[k][nre['file'] == row['file']] = nrz[k][i]
            
        if len(nre) > 0:
        
            print(f'NIRSPec: {len(nre)} extractions')
        
            keys = ['{program} {grating} {filter}'.format(**row) for row in nre]
            un = utils.Unique(keys)
        
            #nre['comment'] = ['targetid={targetid} <br> {spectype} z={z:.4f} &plusmn; {zerr:.4f}'.format(**row) for row in edr]
            nre['tooltip'] = ['{file} z={z:.4f} grade={grade} </br> <img src="https://s3.amazonaws.com/msaexp-nirspec/extractions/{root}/{file}.fnu.png" height=300px/>'.format(**row).replace('.spec.fits','')
                             for row in nre]
            nre['popup'] = ['<a href="https://s3.amazonaws.com/msaexp-nirspec/extractions/{root}/{file}"/> {file} </a> z={z:.4f} grade={grade} </br> <img src="https://s3.amazonaws.com/msaexp-nirspec/extractions/{root}/{file}.fnu.png" height=300px/> <br> <img src="https://s3.amazonaws.com/msaexp-nirspec/extractions/{root}/{file}.flam.png" height=300px/>'.format(**row).replace('.spec.fits.f','.f')
                             for row in nre]
        
            nre['x'], nre['y'] = wcs.all_world2pix(nre['ra'], nre['dec'], 0)
        
            rows.append("var spec_tt = {direction:'auto'};")
            rows.append("var spec_m3 = {color:'#7AE27A',radius:8,weight:2,opacity:0.95,fill:false};")
            rows.append("var spec_m2 = {color:'#E2DF7A',radius:8,weight:2,opacity:0.95,fill:false};")
            rows.append("var spec_m0 = {color:'salmon',radius:8,weight:2,opacity:0.95,fill:false};")
        
            for j, k in enumerate(un.values):
                #kl = 'spec_'+k.lower().replace(' ','_')
                kl = f'nre_{j}'
                rows.append(f'var {kl} = []; // {k}')
                for row in nre[un[k]]:
                
                    if row['grade'] == 3:
                        gm = 'spec_m3'
                    elif row['grade'] == 2:
                        gm = 'spec_m2'
                    else:
                        gm = 'spec_m0'
                    
                    marker = f"""L.circleMarker([{row['y']:.1f},{row['x']:.1f}], {gm}).bindTooltip('{row['tooltip']}', spec_tt).bindPopup('{row['popup']}')"""
                
                    rows.append(f'{kl}.push({marker});')
                            
                rows.append(f"overlays['Spectra {k}'] = L.layerGroup({kl});")
        
    with open(output_file,'w') as fp:
        for row in rows:
            fp.write(row+'\n')
    
    if upload:
        os.system(f'aws s3 cp {output_file} s3://{S3_MAP_PREFIX}/{field}/ --acl public-read')
            