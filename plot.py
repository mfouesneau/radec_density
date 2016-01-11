#!/usr/bin/env python
from __future__ import print_function
import fileinput

# need to parse options such that pylab uses Agg in quiet mode
if __name__ == "__main__":

    opts = (
        ('--title', dict(dest="title", default=None, type='str',
                         help="Figure title")),
        ('-o', '--output', dict(dest="figout", default=None, type='str',
                                help="export figures using this suffixe")),
        ("-q", '--quiet', dict(action="store_false", dest="show",
                               default=True, help="set to skip showing the figures")),
        ("--projection", dict(dest="projection", default='hammer', type='str',
                              help="MAP Coordinate projection base")),
    )

    from optparse import OptionParser
    parser = OptionParser()
    for ko in opts:
        parser.add_option(*ko[:-1], **ko[-1])

    (options, args) = parser.parse_args()

    output = options.__dict__.pop('figout')
    title = options.__dict__.pop('title')
    projection = options.__dict__.pop('projection')

    if len(args) >= 1:
        infile = args[0]
        if infile == '-':
            infile = fileinput.FileInput(args[0])
    else:
        infile = fileinput.FileInput('-')

    if not options.show:
        from matplotlib import use
        use('Agg')

import pylab as plt
from mpl_toolkits.basemap import Basemap
from matplotlib import ticker
import matplotlib.patheffects
import numpy as np
from matplotlib import rcParams


def declare_parula():
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib import cm

    cm_data = [
        [0.2081, 0.1663, 0.5292],
        [0.2116238095, 0.1897809524, 0.5776761905],
        [0.212252381, 0.2137714286, 0.6269714286],
        [0.2081, 0.2386, 0.6770857143],
        [0.1959047619, 0.2644571429, 0.7279],
        [0.1707285714, 0.2919380952, 0.779247619],
        [0.1252714286, 0.3242428571, 0.8302714286],
        [0.0591333333, 0.3598333333, 0.8683333333],
        [0.0116952381, 0.3875095238, 0.8819571429],
        [0.0059571429, 0.4086142857, 0.8828428571],
        [0.0165142857, 0.4266, 0.8786333333],
        [0.032852381, 0.4430428571, 0.8719571429],
        [0.0498142857, 0.4585714286, 0.8640571429],
        [0.0629333333, 0.4736904762, 0.8554380952],
        [0.0722666667, 0.4886666667, 0.8467],
        [0.0779428571, 0.5039857143, 0.8383714286],
        [0.079347619, 0.5200238095, 0.8311809524],
        [0.0749428571, 0.5375428571, 0.8262714286],
        [0.0640571429, 0.5569857143, 0.8239571429],
        [0.0487714286, 0.5772238095, 0.8228285714],
        [0.0343428571, 0.5965809524, 0.819852381],
        [0.0265, 0.6137, 0.8135],
        [0.0238904762, 0.6286619048, 0.8037619048],
        [0.0230904762, 0.6417857143, 0.7912666667],
        [0.0227714286, 0.6534857143, 0.7767571429],
        [0.0266619048, 0.6641952381, 0.7607190476],
        [0.0383714286, 0.6742714286, 0.743552381],
        [0.0589714286, 0.6837571429, 0.7253857143],
        [0.0843, 0.6928333333, 0.7061666667],
        [0.1132952381, 0.7015, 0.6858571429],
        [0.1452714286, 0.7097571429, 0.6646285714],
        [0.1801333333, 0.7176571429, 0.6424333333],
        [0.2178285714, 0.7250428571, 0.6192619048],
        [0.2586428571, 0.7317142857, 0.5954285714],
        [0.3021714286, 0.7376047619, 0.5711857143],
        [0.3481666667, 0.7424333333, 0.5472666667],
        [0.3952571429, 0.7459, 0.5244428571],
        [0.4420095238, 0.7480809524, 0.5033142857],
        [0.4871238095, 0.7490619048, 0.4839761905],
        [0.5300285714, 0.7491142857, 0.4661142857],
        [0.5708571429, 0.7485190476, 0.4493904762],
        [0.609852381, 0.7473142857, 0.4336857143],
        [0.6473, 0.7456, 0.4188],
        [0.6834190476, 0.7434761905, 0.4044333333],
        [0.7184095238, 0.7411333333, 0.3904761905],
        [0.7524857143, 0.7384, 0.3768142857],
        [0.7858428571, 0.7355666667, 0.3632714286],
        [0.8185047619, 0.7327333333, 0.3497904762],
        [0.8506571429, 0.7299, 0.3360285714],
        [0.8824333333, 0.7274333333, 0.3217],
        [0.9139333333, 0.7257857143, 0.3062761905],
        [0.9449571429, 0.7261142857, 0.2886428571],
        [0.9738952381, 0.7313952381, 0.266647619],
        [0.9937714286, 0.7454571429, 0.240347619],
        [0.9990428571, 0.7653142857, 0.2164142857],
        [0.9955333333, 0.7860571429, 0.196652381],
        [0.988, 0.8066, 0.1793666667],
        [0.9788571429, 0.8271428571, 0.1633142857],
        [0.9697, 0.8481380952, 0.147452381],
        [0.9625857143, 0.8705142857, 0.1309],
        [0.9588714286, 0.8949, 0.1132428571],
        [0.9598238095, 0.9218333333, 0.0948380952],
        [0.9661, 0.9514428571, 0.0755333333],
        [0.9763, 0.9831, 0.0538]]

    parula_map = LinearSegmentedColormap.from_list('parula', cm_data)
    cm.register_cmap('parula', cmap=parula_map)
    cm.__dict__['parula'] = cm.get_cmap('parula')
    parula_r_map = LinearSegmentedColormap.from_list('parula_r', cm_data[::-1])
    cm.register_cmap('parula_r', cmap=parula_r_map)
    cm.__dict__['parula_r'] = cm.get_cmap('parula_r')

    parula_map = LinearSegmentedColormap.from_list('parulaW', [[1., 1., 1.]] + cm_data + [[1., 1., 1.]])
    cm.register_cmap('parulaW', cmap=parula_map)
    cm.__dict__['parulaW'] = cm.get_cmap('parulaW')
    parula_r_map = LinearSegmentedColormap.from_list('parulaW_r', ([[1., 1., 1.]] + cm_data + [[1., 1., 1.]])[::-1])
    cm.register_cmap('parulaW_r', cmap=parula_r_map)
    cm.__dict__['parulaW_r'] = cm.get_cmap('parulaW_r')


# configuration
declare_parula()
rcParams['font.family'] = 'serif'
rcParams['text.latex.preamble'] = [r'\usepackage{amsmath}']
rcParams['font.size'] = 14


def get_radec_basemap(projection='moll', resolution='c', lon_0=0, **kwargs):
    """ Generate a Basemap for angular plot

    Parameters
    ----------
    projection: string
        projection name (see: :class:`basemap.Basemap`)

    resolution: string
        resolution of the projection (see: :class:`basemap.Basemap`)

    lon_0: float
        origin of the reference (see: :class:`basemap.Basemap`)

    kwars: dict
        forwarded to :func:`Basemap`

    returns
    -------
    m: Basemap instance
        result from :func:`Basemap`
    """
    # generate the basemap for the plot
    m = Basemap(projection=projection,
                lon_0=lon_0,
                resolution=resolution, **kwargs)

    m.drawmapboundary(fill_color='white')

    m.drawmeridians(np.arange(-180, 180, 360. / 24. * 3.),
                    labels=[0, 0, 0, 0],
                    color='black',
                    dashes=[1,1],
                    labelstyle='+/-',
                    linewidth=0.1)
    m.drawparallels(np.arange(-90, 91, 20),
                    labels=[1,0,0,0],
                    color='black',
                    dashes=[1,1],
                    labelstyle='+/-',
                    linewidth=0.1)

    # add manually the meridian labels for this projection
    ef = matplotlib.patheffects.withStroke(foreground="w", linewidth=3)
    for k in np.arange(-180, 181, 360. / 24. * 3.):
        val = int(k / 360. * 24.)
        x, y = m(k, 0)
        plt.text(x, y, '{0:+d}h'.format(val), path_effects=[ef],
                 ha='center', va='center')

    return m


def plot_ra_dec_histogram2d(density, ra_bins, dec_bins, projection='moll',
                            resolution='c', lon_0=0, cmap=plt.cm.parulaW,
                            vmin=None, vmax=None, basemap=None, **kwargs):
    """
    Plot an histogram2d in a projected system of coordinates

    Parameters
    ----------
    density: ndarray, (nx, ny)
        density map

    ra_bins: ndarray, (nx,)
        ra bin edges

    dec_bins: ndarray, (ny, )
        dec bin edges

    projection: string
        projection name (see: :class:`basemap.Basemap`)

    resolution: string
        resolution of the projection (see: :class:`basemap.Basemap`)

    lon_0: float
        origin of the reference (see: :class:`basemap.Basemap`)

    cmap: plt.colormap
        colormap for :func:`plt.pcolormesh`

    vmin: float
        minimum value to show

    vmax: float
        maximum value to show

    kwars: dict
        forwarded to :func:`plt.pcolormesh`

    returns
    -------
    r: [patches]
        result from :func:`plt.pcolormesh`

    """
    # generate the basemap for the plot
    if basemap is None:
        basemap = get_radec_basemap(projection=projection, lon_0=lon_0,
                                    resolution=resolution)

    dec_bins_2d, ra_bins_2d = np.meshgrid(dec_bins, ra_bins)
    # convert the xs and ys to map coordinates
    xs, ys = basemap(ra_bins_2d, dec_bins_2d)
    if vmin is None:
        vmin = np.nanmin(density)
    if vmax is None:
        vmax = np.nanmax(density)
    r = plt.pcolormesh(xs, ys, density, cmap=cmap, vmin=vmin, vmax=vmax,
                       **kwargs)
    return r


def set_colorbar_MaxNLocator(cb, n):
    """ Set the number of Major Locators on a colorbar """
    tick_locator = ticker.MaxNLocator(nbins=n)
    cb.locator = tick_locator
    cb.update_ticks()


if __name__ == "__main__":

    density = np.loadtxt(infile).T
    nx, ny = density.shape
    dx = 360. / float(nx)
    dy = 180. / float(ny)
    ra_bins = np.arange(-180, 180 + 0.1 * dx, dx)
    dec_bins = np.arange(-90, 90 + 0.1 * dy, dy)

    # projection = "hammer"

    plt.ioff()
    # plot in mollweide projections
    m = get_radec_basemap(projection)
    # from matplotlib.colors import LogNorm
    # r = plot_ra_dec_histogram2d(density + 0.1, ra_bins, dec_bins, basemap=m,
    #                             norm=LogNorm(vmin=1e-3, vmax=density.max())
    #                             )
    r = plot_ra_dec_histogram2d(density, ra_bins, dec_bins, basemap=m)
    # ra, dec = plot.galactic_plane_RADEC()
    # mra, mdec = m(ra, dec)
    # plt.plot(mra, mdec, 'k-', lw=2)
    # add colorbar
    cb = plt.colorbar(r, orientation='horizontal', shrink=0.6)
    cb.set_label('number counts')
    set_colorbar_MaxNLocator(cb, 5)

    if title is not None:
        plt.title(title)
    if output is not None:
        plt.savefig(output)

    if options.show:
        plt.show()
