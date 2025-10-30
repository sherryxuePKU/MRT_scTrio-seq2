#!/bin/python
import loompy as lp
import numpy as np
import scanpy as sc

metadata_path="/path/to/your/sample.csv"
loom_path="/path/to/your/MRT_4000.loom"

x=sc.read_csv(metadata_path)

row_attrs = {
  "Gene": np.array(x.var_names),
}
col_attrs = {
  "CellID": np.array(x.obs_names)
}
lp.create(loom_path,x.X.transpose(),row_attrs,col_attrs)