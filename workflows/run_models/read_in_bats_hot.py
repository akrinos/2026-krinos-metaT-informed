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
from scipy.interpolate import interp1d


bats_data_p = pd.read_csv("../../data/bats/BATS_Bottle_ll_phosphorus_data.csv")
bats_data_p["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in bats_data_p["time"]]
# convert to uM using density of seawater (original units nmol/kg)
bats_data_p["phosphorus"] = [curr/1000*1.025 for curr in bats_data_p["low_level_phosphorus"]]
bats_data_p["site"] = "BATS"

bats_data_n = pd.read_csv("../../data/bats/BATS_Bottle_no3_no2_data.csv")
bats_data_n["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in bats_data_n["time"]]
# convert to uM using density of seawater (original units umol/kg)
bats_data_n["nitrogen"] = [curr*1.025 for curr in bats_data_n["nitrate_nitrite"]]
bats_data_n["site"] = "BATS"

bats_data_temp = pd.read_csv("../../data/bats/BATS_Bottle_temp_data.csv")
bats_data_temp["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in bats_data_temp["time"]]
bats_data_temp["site"] = "BATS"


hot_data_p = pd.read_csv("../../data/hot/HOT_Bottle_ALOHA_ll_phosphorus_data.csv")
hot_data_p["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in hot_data_p["time"]]
# convert to uM using density of seawater (original units nmol/kg)
hot_data_p["phosphorus"] = [curr/1000*1.025 for curr in hot_data_p["llp"]]
hot_data_p["site"] = "HOT"

hot_data_n = pd.read_csv("../../data/hot/HOT_Bottle_ALOHA_ll_nitrogen_data.csv")
hot_data_n["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in hot_data_n["time"]]
# convert to uM using density of seawater (original units nmol/kg)
hot_data_n["nitrogen"] = [curr/1000*1.025 for curr in hot_data_n["lln"]]
hot_data_n["site"] = "HOT"

hot_data_temp = pd.read_csv("../../data/hot/HOT_Bottle_ALOHA_temp_data.csv")
hot_data_temp["DateTime"] = [datetime.fromisoformat(date_curr.replace("Z", "+00:00")) for \
                          date_curr in hot_data_temp["time"]]
hot_data_temp["site"] = "HOT"


temp_data_all = pd.concat([hot_data_temp,bats_data_temp],axis="rows")
nutrient_data_all = pd.concat([bats_data_p.merge(bats_data_n),
                              hot_data_p.merge(hot_data_n)],axis="rows")

nutrient_data_all = nutrient_data_all.loc[nutrient_data_all.depth<5,:].\
    groupby(["DateTime","site"])[["phosphorus","nitrogen"]].mean().\
    reset_index()
