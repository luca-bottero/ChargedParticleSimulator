import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

df = pd.read_csv("output.dat", sep=" ")
df.columns = ["E_0", "Accuracy", "n_timesteps"]

for i,j,k in zip(df["E_0"].to_numpy(), df["Accuracy"].to_numpy(), df["n_timesteps"].to_numpy()):
    print(np.round(i,3), "&", np.round(j, 9), "&", int(k), "\\\\")
