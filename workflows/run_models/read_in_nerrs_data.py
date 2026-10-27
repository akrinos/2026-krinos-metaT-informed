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

nerrs_data = pd.read_csv("../../data/nerrs_data_cleaned.csv")
nerrs_data["DateTime"] = [datetime.strptime(date_curr,"%Y-%m-%d") for \
                          date_curr in nerrs_data["DateTime"]]

nerrs_data_temp = pd.read_csv("../../data/nerrs_data_cleaned_temp.csv")
nerrs_data_temp["DateTime"] = [datetime.strptime(date_curr,"%Y-%m-%d") for \
                          date_curr in nerrs_data_temp["DateTime"]]