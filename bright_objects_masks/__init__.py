from .version import __version__
import lsst.sphgeom as sphgeom
import numpy as np
import healpy as hp
import shapely

__all__ = [
    "call_dc2",
    "generate_masks",
    "radius_study",
    "_multiprocessing.s_multiprocessing",
]



def get_htm_shapelies(depth, return_id=False) :
    """get shapely polygons corresponding to an HTM of given depth
    Parameters
    ----------
    depth : int
        depth of the HTM mapping controlling the resolution of the pixelization

    Returns
    -------
    numpy array
        array of shapely polygons
    list
        IDs of the HTM wedges
    """
    htmMap = sphgeom.HtmPixelization(depth)
    htm_ids = range(*htmMap.universe()[0])

    htm_shapelies = np.array([
        shapely.Polygon(
            np.vstack(
                hp.vec2ang(
                    np.vstack(sphgeom.ConvexPolygon.getVertices(htmMap.triangle(ID))),
                    lonlat=True
                    )
                ).T
            )
        for ID in htm_ids])

    if not return_id :
        return htm_shapelies
    elif return_id :
        return htm_shapelies, htm_ids


def hpix2htm_id(nside, depth, hpixID, buffer=1., nest=True) :
    """get the HTM IDs corresponding to a HEALPix ID
    Parameters
    ----------
    nside : int
        the nside of the HEALPix map
    depth : int
        the depth of the HTM map
    hpixID : int
        HEALPix ID
    buffer : float
        buffer region around HEALPix to consider in arcmin
    nest : bool
        HEALPix ID corresponds to a nested (True) or ringed (False) ordering

    Returns
    -------
    list
        IDs of the HTM wedges overlapping with the HEALPix pixel
    """
    coords = hp.vec2ang(hp.boundaries(nside, hpixID, step=1, nest=nest).T, lonlat=True)
    hpix_poly = shapely.Polygon(np.array(coords).T)
    hpix_poly_buff = hpix_poly.buffer(buffer / 60)

    htmMap = sphgeom.HtmPixelization(depth)
    htm_ids_all = range(*htmMap.universe()[0])
    intersections = np.array([htm_shapely.intersects(hpix_poly_buff) for htm_shapely in get_htm_shapelies(depth)])
    htm_ids_overlap = np.array(htm_ids_all)[intersections]
    
    return htm_ids_overlap


def tract2htm_id(depth, tract_id, skymap, buffer=1.) :
    """get the HTM IDs corresponding to a tract
    Parameters
    ----------
    depth : int
        the depth of the HTM map
    tract_id : int
        ID of the tract
    skymap :
        
    buffer : float
        buffer region around HEALPix to consider in arcmin

    Returns
    -------
    list
        IDs of the HTM wedges overlapping with the HEALPix pixel
    """
    tract_info = skymap[tract_id]
    polygon = tract_info.getOuterSkyPolygon()
    coords = polygon.getVertices()
    tract_poly = shapely.Polygon(np.vstack(hp.vec2ang(np.vstack(coords), lonlat=True)).T)
    tract_poly_buff = tract_poly.buffer(buffer / 60)

    htmMap = sphgeom.HtmPixelization(depth)
    htm_ids_all = range(*htmMap.universe()[0])
    intersections = np.array([htm_shapely.intersects(tract_poly_buff) for htm_shapely in get_htm_shapelies(depth)])
    htm_ids_overlap = np.array(htm_ids_all)[intersections]

    return htm_ids_overlap

