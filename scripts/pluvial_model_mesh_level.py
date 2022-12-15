"""
Module to run the pluvial model at the unstructured mesh level.
"""
import argparse
import json
import os
import time

import dask.array as da
import numpy as np
import pandas as pd
import xarray as xr
import zarr

from scripts.urban.align_grid_size import align_grid_size_improved_urban
from scripts.urban.run_scs_mesh_level import run_scs
from scripts.util.connectionsConfig import (
    GCS_BUCKET,
    GCS_URBAN_STATIC_PATH,
    INPUT_DIR,
    OUTPUT_DIR,
    TEMP_DIR,
)
from scripts.util.get_logger import getLogger

logger = getLogger(__name__)


def run_urban_model_mesh_level(run_id: str, urban_static_data_version: str) -> None:
    """
    Run urban model at the mesh level
    """
    # Read in the Urban Pipeline run arguments
    urban_pipeline_runs_args_file = os.path.join(
        INPUT_DIR, run_id, "urban", "urban_pipeline_run_args.json"
    )
    with open(urban_pipeline_runs_args_file, "r", encoding="utf-8") as f:
        run_args = json.load(f)
    hot_start_precip_path = run_args["hot_start_precip"]
    hot_start_infil_path = run_args["hot_start_infil"]
    scs_tmp_path = os.path.join(TEMP_DIR, run_id, "urban")
    scs_output_path = os.path.join(OUTPUT_DIR, run_id, "urban")

    t1 = time.time()
    urban_map_file = os.path.join(
        GCS_BUCKET,
        GCS_URBAN_STATIC_PATH,
        f"urban_grid_map_{urban_static_data_version}.zarr",
    )
    logger.info(f"Reading urban-schism mapping file: {urban_map_file}")
    ds_hires = xr.open_zarr(
        urban_map_file, consolidated=True, storage_options={"anon": True}
    )
    # Limit the map to only hires cells (with city boundaries)
    hires_mask = ds_hires.resolution.values == 1
    initial_abstraction_values = ds_hires.initial_abstraction.values[hires_mask]
    cell_area_values = ds_hires.area.values[hires_mask]
    logger.info(f"Read in ds_hires in {(time.time() - t1):.2f} secs")

    # This gets created in preprocess_urban_rainfall.py
    df_forecast = pd.read_csv(f"{scs_tmp_path}/forecast_data_sources_{run_id}.csv")
    logger.info("Data sources and paths:")
    for tpl in df_forecast.itertuples():
        logger.info(f"{tpl.data_source}; {tpl.filepath}")

    infil_dir = os.path.join(scs_tmp_path, "infil")
    runoff_dir = os.path.join(scs_output_path, "runoff")
    converter_output_dir = os.path.join(scs_output_path, "urban_to_schism")

    if not os.path.exists(infil_dir):
        os.mkdir(infil_dir)
    if not os.path.exists(runoff_dir):
        os.mkdir(runoff_dir)
    if not os.path.exists(converter_output_dir):
        os.mkdir(converter_output_dir)

    t1 = time.time()
    wb_mask = ds_hires.waterbody.values[hires_mask] == 1
    hires_x_vals = ds_hires.x.values[hires_mask]
    hires_y_vals = ds_hires.y.values[hires_mask]

    # First do the regridding of rainfall to high-res mesh
    gps = df_forecast.groupby(by="data_source")
    lst_regridded_rainfall_files = []
    for _, df in gps:
        for i, tpl in enumerate(df.itertuples()):
            logger.info(
                f"Regridding {os.path.basename(tpl.filepath)}; data_source: {tpl.data_source}"
            )
            source_xyz = zarr.load(tpl.filepath)
            source_coords = source_xyz[:, 0:2]
            if i == 0:  # only need to do this once per data source
                t1 = time.time()
                nn_inds = align_grid_size_improved_urban(
                    source_coords, hires_x_vals, hires_y_vals
                )
                regridded_rainfall_array = source_xyz[:, 2][nn_inds]
                logger.info(
                    f"Regridded with align_grid_size in: {(time.time() - t1):.2f} secs"
                )
            else:
                regridded_rainfall_array = source_xyz[:, 2][nn_inds]
            # Save the regridded file
            regridded_path = f"{tpl.filepath[:-4]}_highres.zarr"
            da_arr = da.from_array(regridded_rainfall_array, chunks=1000000)
            da_arr.to_zarr(regridded_path, overwrite=True)
            lst_regridded_rainfall_files.append(regridded_path)
    logger.info(
        f"Regridded the rainfall for this forecast in: {(time.time() - t1) / 60:.2f} mins"
    )
    # Make sure they're ordered by time
    lst_regridded_rainfall_files.sort()

    # Assemble infiltration and runoff filepaths
    lst_infil_paths = []
    lst_runoff_paths = []
    for tpl in df_forecast.itertuples():
        filename = os.path.basename(tpl.filepath)
        infil_path = f"{infil_dir}/{filename[:-4]}_infil.zarr"
        lst_infil_paths.append(infil_path)
        runoff_path = f"{runoff_dir}/{filename[:-4]}_runoff.zarr"
        lst_runoff_paths.append(runoff_path)

    # Get hostarts if available
    t1 = time.time()
    if (hot_start_precip_path is None) | (hot_start_infil_path is None):
        precip_previous = np.zeros(initial_abstraction_values.shape)
        infil_previous = np.zeros(initial_abstraction_values.shape)
        logger.info("Precip and infil previous are set to zero")
    else:
        precip_previous = da.from_zarr(hot_start_precip_path).compute()
        infil_previous = da.from_zarr(hot_start_infil_path).compute()
        logger.info("Read in rainfall and infil hotstarts")

    # Now loop over each step sequentially and calculate runoff
    num_timesteps = len(df_forecast.index)
    for i in range(num_timesteps):
        logger.info(f"Calculating runoff (cms) for: {lst_regridded_rainfall_files[i]}")
        target_precip = zarr.load(lst_regridded_rainfall_files[i])
        precip, infil = run_scs(
            initial_abstraction_values,
            target_precip,
            wb_mask,
            precip_previous,
            infil_previous,
            cell_area_values,
            lst_runoff_paths,
            i,
        )
        precip_previous = precip
        infil_previous = infil
        da_arr = da.from_array(infil, chunks=1000000)
        da_arr.to_zarr(lst_infil_paths[i], overwrite=True)
    logger.info(
        f"Completed CN method runoff calcs and conversion to cms in: {(time.time() - t1) / 60:.2f} mins"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_id")
    parser.add_argument("urban_static_data_version")

    args = parser.parse_args()
    logger.info("Args: %s", args)

    run_urban_model_mesh_level(args.run_id, args.urban_static_data_version)
