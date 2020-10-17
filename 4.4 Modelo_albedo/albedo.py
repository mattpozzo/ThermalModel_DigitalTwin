############## Albedo Model Block ##########################
  # Coded by Universitwin - Digital Twin Team
  # April 2020
  # Started by: Nicolás Valentín Conde
  # Version: V02
  # Latest version by:
  # Latest version date:

# Libraries...
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Functions

# Input

# Input 1: Albedo map
df = pd.read_csv('../4.8 CSVs/albedo.CSV')
df = df.iloc[:,1:]
df[df == df.where(df > 1)] = 0

# Input 2: Position, attitude, time

# Input 3: Experiment configuration

# Algorithm...
# Assign to each face of experiment, a value of albedo
# for every position.

# Output...

# Visualization
  # Input 1 visualization
d1 = np.arange(0, 1799, 1)
d2 = np.arange(0, 3599, 1)
X, Y = np.meshgrid(d1, d2)
plt.subplots()
plt.imshow(df, cmap='bone')
plt.colorbar()
plt.show()

########## ##     ####     ######
    ##     ##    ##  ##    ##   ##
    ##     ##   ##    ##   ##   ##
    ##     ##   ########   ######
    ##     ##  ##      ##  ##   ## 
    ##     ##  ##      ##  ##   ##
    ##     ##  ##      ##  ######
