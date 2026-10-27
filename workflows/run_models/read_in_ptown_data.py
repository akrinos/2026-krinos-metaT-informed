import os
import pandas as pd
import numpy as np
import scipy
import math
import plotnine
import matplotlib.pyplot as plt
import warnings
from scipy.integrate import solve_ivp
from datetime import datetime
from datetime import timedelta

from model_functions import *

ammon_ptown_1 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2000_2021.csv") 
ammon_ptown_1["Date"] = [datetime.strptime(date_string, "%d-%b-%Y %H:%M:%S") for date_string in ammon_ptown_1.STAT_ARRIV]
ammon_ptown_2 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2022_2023.csv")
ammon_ptown_2["Date"] = [datetime.strptime(date_string, "%Y-%m-%d") for date_string in ammon_ptown_2.STAT_ARRIV]

ammon_ptown = pd.concat([ammon_ptown_1,ammon_ptown_2])

ammon_ptown = ammon_ptown.loc[[param_code in ["NH4"]\
                             for param_code in ammon_ptown["PARAM_CODE"]],:].\
    groupby(["PARAM_CODE","STAT_ARRIV","STAT_ID","Date"]).FFU_RESULT.mean().reset_index().\
    rename({"FFU_RESULT":"ammonium"},axis="columns")

ammon_ptown = ammon_ptown.groupby(["Date","STAT_ID"]).ammonium.mean().reset_index()

ammon_ptown["DateMonthYr"] = [dt.replace(day=1) for dt in ammon_ptown["Date"]]
ammon_ptown = ammon_ptown.sort_values("DateMonthYr")
ammon_ptown["Month"] = [dt.month for dt in ammon_ptown["Date"]]
ammon_ptown["Year"] = [dt.year for dt in ammon_ptown["Date"]]

time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-11-01 00:00","%Y-%m-%d %H:%M"), freq='M')
to_infill = pd.DataFrame({"DateTime":time_vect})
to_infill["Year"] = [dt.year for dt in to_infill.DateTime]
to_infill["Month"] = [dt.month for dt in to_infill.DateTime]
to_infill["Hour"] = [dt.hour for dt in to_infill.DateTime]
to_infill = to_infill.merge(ammon_ptown,on=["Month","Year"],how="left")
ammon_1 = infill_data(to_infill,"ammonium")
time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-10-30 00:00","%Y-%m-%d %H:%M"), freq='D')
time_frame_2 = pd.DataFrame({"DateTime":time_vect})
time_frame_2["Year"] = [dt.year for dt in time_frame_2.DateTime]
time_frame_2["Month"] = [dt.month for dt in time_frame_2.DateTime]
time_frame_2["Hour"] = [dt.hour for dt in time_frame_2.DateTime]
ammon = time_frame_2.merge(ammon_1.loc[:,["Month","Year","ammonium"]].drop_duplicates(subset=["Month","Year"]),
                     on=["Month","Year"],how="left")

nitr_ptown_1 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2000_2021.csv") 
nitr_ptown_1["Date"] = [datetime.strptime(date_string, "%d-%b-%Y %H:%M:%S") for date_string in nitr_ptown_1.STAT_ARRIV]
nitr_ptown_2 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2022_2023.csv")
nitr_ptown_2["Date"] = [datetime.strptime(date_string, "%Y-%m-%d") for date_string in nitr_ptown_2.STAT_ARRIV]

nitr_ptown = pd.concat([nitr_ptown_1,nitr_ptown_2])
nitr_ptown = nitr_ptown.loc[nitr_ptown["DEPTH"]<5,:]

nitr_ptown = nitr_ptown.loc[[param_code in ["NO2","NO3","NO3+NO2"]\
                             for param_code in nitr_ptown["PARAM_CODE"]],:].\
    groupby(["PARAM_CODE","STAT_ARRIV","STAT_ID","Date"]).FFU_RESULT.mean().reset_index()

nitr_ptown = nitr_ptown.pivot(index=["STAT_ARRIV","STAT_ID","Date"],
                              columns="PARAM_CODE", values="FFU_RESULT").\
    rename({"NO3+NO2":"NO3_NO2"},axis="columns").reset_index()
nitr_ptown["nitrogen"] = [currno3_no2 if currno3_no2==currno3_no2 else\
                          currno3+currno2 for currno3_no2,currno3,currno2 in\
                          zip(nitr_ptown["NO3_NO2"],nitr_ptown["NO3"],
                              nitr_ptown["NO2"])]
nitr_ptown = nitr_ptown.groupby(["Date","STAT_ID"]).nitrogen.mean().reset_index()

nitr_ptown["DateMonthYr"] = [dt.replace(day=1) for dt in nitr_ptown["Date"]]
nitr_ptown = nitr_ptown.sort_values("DateMonthYr")
nitr_ptown["Month"] = [dt.month for dt in nitr_ptown["Date"]]
nitr_ptown["Year"] = [dt.year for dt in nitr_ptown["Date"]]

time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-11-01 00:00","%Y-%m-%d %H:%M"), freq='M')
to_infill = pd.DataFrame({"DateTime":time_vect})
to_infill["Year"] = [dt.year for dt in to_infill.DateTime]
to_infill["Month"] = [dt.month for dt in to_infill.DateTime]
to_infill["Hour"] = [dt.hour for dt in to_infill.DateTime]
to_infill = to_infill.merge(nitr_ptown,on=["Month","Year"],how="left")
nitr_1 = infill_data(to_infill,"nitrogen")
time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-10-30 00:00","%Y-%m-%d %H:%M"), freq='D')
time_frame_2 = pd.DataFrame({"DateTime":time_vect})
time_frame_2["Year"] = [dt.year for dt in time_frame_2.DateTime]
time_frame_2["Month"] = [dt.month for dt in time_frame_2.DateTime]
time_frame_2["Hour"] = [dt.hour for dt in time_frame_2.DateTime]
nitr = time_frame_2.merge(nitr_1.loc[:,["Month","Year","nitrogen"]].drop_duplicates(subset=["Month","Year"]),
                     on=["Month","Year"],how="left")
nitr = nitr.merge(ammon)
nitr["nitrogen"]= [curr_nit + curr_amm for curr_nit,curr_amm in zip(nitr["nitrogen"],nitr["ammonium"])]

phos_ptown_1 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2000_2021.csv") 
phos_ptown_1["Date"] = [datetime.strptime(date_string, "%d-%b-%Y %H:%M:%S") for date_string in phos_ptown_1.STAT_ARRIV]
phos_ptown_2 = pd.read_csv("../../data/nutrients/CC_wet_chemistry_2022_2023.csv")
phos_ptown_2["Date"] = [datetime.strptime(date_string, "%Y-%m-%d") for date_string in phos_ptown_2.STAT_ARRIV]

phos_ptown = pd.concat([phos_ptown_1,phos_ptown_2])
phos_ptown = phos_ptown.loc[phos_ptown["PARAM_CODE"]=="PO4",:]
phos_ptown = phos_ptown.loc[phos_ptown["DEPTH"]<5,:]
    
phos_ptown = phos_ptown.groupby(["Date","STAT_ID"]).FFU_RESULT.mean().reset_index().\
    rename({"FFU_RESULT":"phosphorus"},axis="columns")

phos_ptown["DateMonthYr"] = [dt.replace(day=1) for dt in phos_ptown["Date"]]
phos_ptown = phos_ptown.sort_values("DateMonthYr")
phos_ptown["Month"] = [dt.month for dt in phos_ptown["Date"]]
phos_ptown["Year"] = [dt.year for dt in phos_ptown["Date"]]

time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-11-01 00:00","%Y-%m-%d %H:%M"), freq='M')
to_infill = pd.DataFrame({"DateTime":time_vect})
to_infill["Year"] = [dt.year for dt in to_infill.DateTime]
to_infill["Month"] = [dt.month for dt in to_infill.DateTime]
to_infill["Hour"] = [dt.hour for dt in to_infill.DateTime]
to_infill = to_infill.merge(phos_ptown,on=["Month","Year"],how="left")
phos_1 = infill_data(to_infill,"phosphorus")
time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-10-30 00:00","%Y-%m-%d %H:%M"), freq='D')
time_frame_2 = pd.DataFrame({"DateTime":time_vect})
time_frame_2["Year"] = [dt.year for dt in time_frame_2.DateTime]
time_frame_2["Month"] = [dt.month for dt in time_frame_2.DateTime]
time_frame_2["Hour"] = [dt.hour for dt in time_frame_2.DateTime]
phos = time_frame_2.merge(phos_1.loc[:,["Month","Year","phosphorus"]].drop_duplicates(subset=["Month","Year"]),
                     on=["Month","Year"],how="left")
                                     
sili_ptown = pd.read_csv("../../data/silica.csv")

sili_ptown["Date"] = [datetime.strptime(date_string, "%d-%b-%Y %H:%M:%S") for date_string in sili_ptown.STAT_ARRIV]
sili_ptown = sili_ptown.groupby(["Date","STAT_ID"]).FFU_RESULT.mean().reset_index().\
    rename({"FFU_RESULT":"silica"},axis="columns")

sili_ptown["DateMonthYr"] = [dt.replace(day=1) for dt in sili_ptown["Date"]]
sili_ptown = sili_ptown.sort_values("DateMonthYr")
sili_ptown["Month"] = [dt.month for dt in sili_ptown["Date"]]
sili_ptown["Year"] = [dt.year for dt in sili_ptown["Date"]]

time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-11-01 00:00","%Y-%m-%d %H:%M"), freq='M')
to_infill = pd.DataFrame({"DateTime":time_vect})
to_infill["Year"] = [dt.year for dt in to_infill.DateTime]
to_infill["Month"] = [dt.month for dt in to_infill.DateTime]
to_infill["Hour"] = [dt.hour for dt in to_infill.DateTime]
to_infill = to_infill.merge(sili_ptown,on=["Month","Year"],how="left")
sili_1 = infill_data(to_infill,"silica")
time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-10-30 00:00","%Y-%m-%d %H:%M"), freq='D')
time_frame_2 = pd.DataFrame({"DateTime":time_vect})
time_frame_2["Year"] = [dt.year for dt in time_frame_2.DateTime]
time_frame_2["Month"] = [dt.month for dt in time_frame_2.DateTime]
time_frame_2["Hour"] = [dt.hour for dt in time_frame_2.DateTime]
sili = time_frame_2.merge(sili_1.loc[:,["Month","Year","silica"]].drop_duplicates(subset=["Month","Year"]),
                     on=["Month","Year"],how="left")


### TEMPERATURE CALCULATION ###
time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-11-01 00:00","%Y-%m-%d %H:%M"), freq='M')

temps_ptown = pd.read_csv("../../data/temperatures.csv")

temps_ptown["Date"]=[datetime.strptime(date_string, "%Y-%m-%d") for date_string in temps_ptown.Date]
temps_ptown = temps_ptown.groupby(["Date","STAT_ID"]).FFU_RESULT.mean().reset_index().\
    rename({"FFU_RESULT":"temperature"},axis="columns")

temps_ptown["DateMonthYr"] = [dt.replace(day=1) for dt in temps_ptown["Date"]]
temps_ptown = temps_ptown.sort_values("DateMonthYr")
temps_ptown["Month"] = [dt.month for dt in temps_ptown["Date"]]
temps_ptown["Year"] = [dt.year for dt in temps_ptown["Date"]]

to_infill = pd.DataFrame({"DateTime":time_vect})
to_infill["Year"] = [dt.year for dt in to_infill.DateTime]
to_infill["Month"] = [dt.month for dt in to_infill.DateTime]
to_infill["Hour"] = [dt.hour for dt in to_infill.DateTime]
to_infill = to_infill.merge(temps_ptown,on=["Month","Year"],how="left")
temp_1 = infill_data(to_infill,"temperature")

time_vect = pd.date_range(datetime.strptime('2000-02-01 00:00',"%Y-%m-%d %H:%M"),
                          datetime.strptime("2023-10-30 00:00","%Y-%m-%d %H:%M"), freq='h')

time_frame_2 = pd.DataFrame({"DateTime":time_vect})
time_frame_2["Year"] = [dt.year for dt in time_frame_2.DateTime]
time_frame_2["Month"] = [dt.month for dt in time_frame_2.DateTime]
time_frame_2["Day"] = [dt.day for dt in time_frame_2.DateTime]
time_frame_2["Hour"] = [dt.hour for dt in time_frame_2.DateTime]
daily_temp_fr_2 = time_frame_2.merge(temp_1.loc[:,["Month","Year","temperature"]].drop_duplicates(subset=["Month","Year"]),
                     on=["Month","Year"],how="left")

daily_temp_fr_2["temperature_hourly"]=[temperature + np.sin(2*np.pi*(Hour-6)/24)/2 for temperature,Hour in\
                                       zip(daily_temp_fr_2.temperature,daily_temp_fr_2.Hour)]
daily_temp_fr_2["temperature_daily"] = daily_temp_fr_2["temperature"]

daily_temp_fr_2 = daily_temp_fr_2.groupby(["Year","Month","Day","DateTime"])[["temperature_hourly","temperature_daily"]].mean().reset_index()
#########