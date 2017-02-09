# -*- coding: utf-8 -*-
# To run from command terminal:
# Start > Run > cmd
# for non-gg species:
#   G: & cd G:\Team_Folders\Steph\scripts & C:\Python27\ArcGIS10.3\python.exe run_sp.py runs_todo_gg_2015.csv
# for gg species, add 2nd argument:
#   G: & cd G:\Team_Folders\Steph\scripts & C:\Python27\ArcGIS10.3\python.exe run_sp.py runs_todo_gg_xxxx.csv gg

# Import arcpy module
import arcgisscripting, shutil, os, csv, sys

# get input arguments
runs_csv = sys.argv[1] # "G:\\Team_Folders\\Steph\\scripts\\runs_todo.csv"
if len(sys.argv) > 2:
    dest_flag = sys.argv[2] # "G:\\Team_Folders\\Steph\\scripts\\runs_todo.csv"
    # assume gg
    dest_ids = "2491;2415;2338;2261;2260;2186;2185;2335;2109;2108;2032;1956;1883;1812;1811;1743;1611;1546;1545;1481;1416;1352;1351;1290;1168"
else:
    dest_ids = "1680;1681;1615;1616"

gp = arcgisscripting.create()
gp.AddToolbox('C:\\Program Files\\GeoEco\\ArcGISToolbox\\Marine Geospatial Ecology Tools.tbx', 'GeoEco')

watermask_tif = "G:\\Data_Layers\\study_area\\ecorgn_water_ga83e_steph_val.tif"
cache_dir     = "G:\\Team_Folders\\MGETdata\\MGET_Cache"

# read in csv with runs per row and variables per column
csv_file = open(runs_csv, 'r')
runs = csv.DictReader(csv_file)

# iterate over runs
for run in runs:
    # assign variables
    for var in run.keys():
        globals()[var] = run[var]
    
    # set species-specific paths
    print sp
    reefid_tif     = "G:\\Team_Folders\\Steph\\nda_{species}\\{species}_reefid.tif".format(species=sp)
    pctcover_tif   = "G:\\Team_Folders\\Steph\\nda_{species}\\{species}_pct_m_n.tif".format(species=sp)
    #simulation_dir = "G:\\Team_Folders\\Steph\\{prefix}_{species}_{suffix}_simulation".format(prefix=pfx, species=sp, suffix=sfx)
    simulation_dir = "G:\\Team_Folders\\Steph\\{species}_{year}\\{prefix}_{species}_{suffix}_simulation".format(year=yr, prefix=pfx, species=sp, suffix=sfx)
    #results_dir    = "G:\\Team_Folders\\Steph\\{prefix}_{species}_{suffix}_results".format(prefix=pfx, species=sp, suffix=sfx)
    results_dir    = "G:\\Team_Folders\\Steph\\{species}_{year}\\{prefix}_{species}_{suffix}_results".format(year=yr, prefix=pfx, species=sp, suffix=sfx)
        
    # Process: Create Larval Dispersal Simulation From ArcGIS Rasters
    print '  LarvalDispersalCreateSimulationFromArcGISRasters_GeoEco ->', simulation_dir
    if os.path.isdir(simulation_dir):
        shutil.rmtree(simulation_dir)
    gp.LarvalDispersalCreateSimulationFromArcGISRasters_GeoEco(
        simulation_dir, reefid_tif, pctcover_tif, watermask_tif, "false")

    # Process: Load HYCOM GLBa0.08 Currents Into Larval Dispersal Simulation
    print '  LarvalDispersalLoadHYCOMGLBa0084DEquatorialCurrentsIntoSimulation_GeoEco ->', cache_dir
    endDate = datetime.datetime.strptime(startDate, "%m/%d/%Y") + datetime.timedelta(days=int(durationDays))
    gp.LarvalDispersalLoadHYCOMGLBa0084DEquatorialCurrentsIntoSimulation_GeoEco(
        simulation_dir, startDate, endDate, "0", "false", "CUBIC", "Del2a", "120", "240", cache_dir)

    # Process: Run Larval Dispersal Simulation (2012 Algorithm)
    print '  LarvalDispersalRunSimulation2012_GeoEco ->', results_dir
    if os.path.isdir(results_dir):
        shutil.rmtree(results_dir)
    os.mkdir(results_dir)
    gp.LarvalDispersalRunSimulation2012_GeoEco(
        simulation_dir, results_dir, startDate, durationDays, "0.5", "24", "2000", "0.01", "0.8", "false", "", dest_ids, "", "50")

    # Process: Visualize Larval Dispersal Simulation Results (2012 Algorithm)
    print '  LarvalDispersalVisualizeResults2012_GeoEco ->', results_dir
    gp.LarvalDispersalVisualizeResults2012_GeoEco(
        simulation_dir, results_dir, "output.gdb", "", "", "true", "0.00001", "false", "true", "0.00001", "Quantity")

csv_file.close()