import os
import numpy as np

def draw_axis_labels(ax=None, nlabel=3, format='latex'):
    """
    Draw rounded axis labels in DMS format
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.ticker import MultipleLocator, FixedLocator
    #from descartes import PolygonPatch
    from sregion import patch_from_polygon
    
    from astropy import units as u
    from astropy.coordinates import angles
    
    if ax is None:
        ax = plt.gca()
    
    dy = np.abs(np.diff(ax.get_ylim()))
    dx = np.abs(np.diff(ax.get_xlim()))
    
    yopts = [30./60, 1,2,4,8,16,30,60, 120, 240]
    yminor = [5./60, 10./60,1,2,2,4,10,15, 30, 60]
    ymaj = np.maximum(1., np.round(dy*60/nlabel))
    y_i = int(np.round(np.interp(ymaj, yopts, range(len(yopts)), left=0)))
    ymaj = yopts[y_i]
    yminor = yminor[y_i]
    
    xopts =  [2, 5, 10, 15,30,60,120, 240, 480]
    xminor = [1, 1, 5, 5, 10, 15, 30, 60, 120]
    xmaj = np.maximum(1, np.round(dx*24/360*3600/nlabel))
    x_i = int(np.round(np.interp(xmaj, xopts, range(len(xopts)), left=0)))
    xmaj = xopts[x_i]
    xminor = xminor[x_i]
    
    ax.xaxis.set_major_locator(MultipleLocator(xmaj/3600*360/24))
    ax.yaxis.set_major_locator(MultipleLocator(ymaj/60))

    ax.xaxis.set_minor_locator(MultipleLocator(xminor/3600*360/24))
    ax.yaxis.set_minor_locator(MultipleLocator(yminor/60))
    
    xcoo = [angles.Longitude(t*u.deg) for t in ax.get_xticks()]
    ycoo = [angles.Latitude(t*u.deg) for t in ax.get_yticks()]
    
    ax.xaxis.set_major_locator(FixedLocator(ax.get_xticks()))
    ax.yaxis.set_major_locator(FixedLocator(ax.get_yticks()))
    
    ax.set_xticklabels([t.to_string(u.hourangle,
                                    pad=True, fields=2+((xmaj < 60)),
                                    precision=0, format=format)
                        for i, t in enumerate(xcoo)])
                        
    ax.set_yticklabels([t.to_string(u.deg, pad=True, fields=2+(ymaj < 1), 
                        format=format) for t in ycoo])


def insert_legacysurveys_thumbnail(ax, pixscale=1, layers=['ls-dr9', 'sdss', 'unwise-neo7'], zorder=-1000, alpha=0.8, aspect='auto', verbose=True, grayscale=False, grayscale_params=[99, 1.5, -0.02], **kwargs):
    """
    Insert thumbnail from LegacySurveys into a plot axis with ra/dec celestial
    coordinates (see https://www.legacysurvey.org/viewer/urls)
    
    Parameters
    ----------
    ax : `~matplotlib.axes._subplots.AxesSubplot`
        Axis to hold the thumbnail
    
    pixscale : int
        Pixel scale
    
    layers : list
        Layer to pull.  The script will cycle through ``layers`` until a valid
        layer is found
    
    zorder : int
        zorder level on ``ax`` to `imshow` the thumbnail
    
    alpha : float
        Transparency alpha
    
    aspect : float, 'auto'
        Aspect ratio for `imshow`.  
    
    grayscale : bool
        If True, plot as reverse grayscale summed from RGB thumbnail.  If > 1
        then don't reverse grayscale
    
    grayscale_params : [max_percentile, scale_max, scale_min]
        Scale grayscale image with
        `vmax=np.percentile(img, max_percentile)*scale_max` and 
        `vmin=scale_min*vmax`
        
    Returns
    -------
    url : str
        Thumbnail URL
    
    img : array
        Thumbnail array
    
    And pulls the thumbnail and inserts it into `ax`
    
    """
    import numpy as np
    
    try:
        from skimage import io
    except ImportError:
        print('insert_legacysurveys_thumbnails: Failed to import skimage')
        return '', None
        
    from urllib.error import HTTPError
    
    ls_url = "https://www.legacysurvey.org/viewer/jpeg-cutout?"
    ls_url += "ra={ra}&dec={dec}&pixscale={ps}&width={sw}&height={sh}"
    ls_url += "&layer={layer}"
    
    xl = ax.get_xlim()
    yl = ax.get_ylim()
                
    dw = np.abs(xl[1]-xl[0])
    sw = int(np.round(dw*np.cos(np.mean(yl)/180*np.pi)*3600/pixscale))
    sh = int(np.round((yl[1]-yl[0])*3600/pixscale))
    
    img = None
    for layer in layers:
        url = ls_url.format(ra=np.mean(xl), dec=np.mean(yl),
                            sw=sw, sh=sh, ps=pixscale, layer=layer)
        try:
            img = io.imread(url)
            if verbose:
                print(url)
            break
        except HTTPError:
            continue

    if img is not None:
        img = np.flipud(img)
        if xl[0] < xl[1]:
            print('insert_legacysurveys_thumbnail: RA axis is flipped')
            img = np.fliplr(img)
        
        if grayscale:
            img = (img*1.).sum(axis=2)
            if grayscale == 1:
                cmap = 'gray_r'
            else:
                cmap = 'gray'
            
            img -= np.median(img)
            vma = np.percentile(img, grayscale_params[0])*grayscale_params[1]
            vmi = grayscale_params[2]*vma
            
            ax.imshow(img, extent=xl+yl, zorder=zorder, alpha=alpha, 
                      aspect='auto', cmap=cmap, vmin=vmi, vmax=vma)
             
        else:
            ax.imshow(img, extent=xl+yl, zorder=zorder, alpha=alpha, 
                      aspect='auto')
    
    return url, img
