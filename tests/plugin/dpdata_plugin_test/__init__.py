from __future__ import annotations

import numpy as np

from dpdata.data_type import Axis, DataType, register_data_type

# test data type

register_data_type(
    DataType("foo", np.ndarray, (Axis.NFRAMES, 2, 4), required=False), labeled=True
)

register_data_type(
    DataType("foo", np.ndarray, (Axis.NFRAMES, 3, 3), required=False), labeled=False
)

register_data_type(
    DataType("bar", np.ndarray, (Axis.NFRAMES, Axis.NATOMS, -1), required=False),
    labeled=True,
)

ep = None
