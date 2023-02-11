#!/usr/bin/env python3
'''
Shared functions for Xarray
'''

#########################
# FILE FORMAT FUNCTIONS #
#########################


def xr_preproc(ds):
    """Preprocessing functuin to adjust coordinate and dimensions
    names to a common format. To be called by xr.open_mf_dataset()"""

    rename_dict = {
        "time_counter": "time",
        "pressure_levels": "plev",
        "plevel": "plev",
        "longitude": "lon",
        "nav_lon": "lon",
        "latitude": "lat",
        "nav_lat": "lat",
        "values": "cell",
    }

    # compact renaming
    for old_name in [name for name in rename_dict if name in ds.dims or name in ds.coords]:
        ds = ds.rename({old_name: rename_dict[old_name]})

    # fix for NEMO eORCA grid (nav_lon, nav_lat)
    for h in ['lon', 'lat']:
        for f in ['', 'grid_T']:
            g = 'nav_' + h + '_' + f
            if g in list(ds.coords):
                ds = ds.rename({g: h})

    # fix for NEMO eORCA grid (x_grid_T, etc.)
    for h in ['x', 'y']:
        for f in ['grid_T']:
            g = h + '_' + f
            if g in list(ds.dims):
                ds = ds.rename({g: h})

    # # print(ds)
    # if 'time_counter' in list(ds.dims):
    #     ds = ds.rename({"time_counter": "time"})

    # if 'time_counter' in list(ds.coords):
    #     ds = ds.rename({"time_counter": "time"})

    # if 'pressure_levels' in list(ds.coords):
    #     ds = ds.rename({"pressure_levels": "plev"})

    # if 'plevel' in list(ds.dims):
    #     ds = ds.rename({"plevel": "plev"})

    # if 'longitude' in list(ds.dims):
    #     ds = ds.rename({"longitude": "lon"})

    # if 'latitude' in list(ds.dims):
    #     ds = ds.rename({"latitude": "lat"})

    # if 'longitude' in list(ds.coords):
    #     ds = ds.rename({"longitude": "lon"})

    # if 'nav_lon' in list(ds.coords):
    #     ds = ds.rename({"nav_lon": "lon"})

    # if 'latitude' in list(ds.coords):
    #     ds = ds.rename({"latitude": "lat"})

    # if 'nav_lat' in list(ds.coords):
    #     ds = ds.rename({"nav_lat": "lat"})

    # if 'values' in list(ds.dims):
    #     ds = ds.rename({"values": "cell"})

    return ds


def adjust_clim_file(cfield, remove_zero=False):
    """Routine to fix file format of climatology"""

    # fix coordinates
    org = ['LONGITUDE', 'LATITUDE', 'lev']
    new = ['lon', 'lat', 'plev']
    for o, n in zip(org, new):
        if o in cfield.coords:
            cfield = cfield.rename({o: n})

    # extract data_array
    cname = list(cfield.data_vars)[-1]
    field = cfield[cname]

    if remove_zero:
        field = field.where(field != 0)

    # convert vertical levels
    if 'plev' in cfield.coords:
        field = field.metpy.convert_coordinate_units('plev', 'Pa')

    return field
