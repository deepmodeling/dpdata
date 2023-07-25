import numpy as np

from dpdata.data_type import Axis, DataType, register_data_type

# test data type

register_data_type(
    DataType("foo", np.ndarray, (Axis.NFRAMES, 2, 4), required=False), labeled=True
)

ep = None
