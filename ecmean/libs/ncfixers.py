#!/usr/bin/env python3
'''
Shared functions for Xarray
'''

#########################
# FILE FORMAT FUNCTIONS #
#########################


def xr_preproc(ds):
    """Preprocessing function to adjust coordinate and dimension
    names to a common format. To be called by xr.open_mf_dataset()

    Parameters:
    ds (xarray.Dataset): The dataset to be preprocessed.

    Returns:
    xarray.Dataset: The preprocessed dataset with adjusted coordinate and dimension names.
    """

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

    # safe check for NEMO output in domain_cfg.nc
    if 'nav_lon' in ds.data_vars and 'nav_lat' in ds.data_vars:
        ds = ds.set_coords(['nav_lon', 'nav_lat'])

    # compact renaming
    for old_name in [name for name in rename_dict if name in ds.dims or name in ds.coords]:
        ds = ds.rename({old_name: rename_dict[old_name]})

    # fix for NEMO eORCA grid (nav_lon, nav_lat)
    for h in ['lon', 'lat']:
        for f in ['', 'grid_T']:
            g = 'nav_' + h + '_' + f
            if g in ds.coords:
                ds = ds.rename({g: h})

    # fix for NEMO eORCA grid (x_grid_T, etc.)
    for h in ['x', 'y']:
        for f in ['grid_T']:
            g = h + '_' + f
            if g in ds.dims:
                ds = ds.rename({g: h})

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
