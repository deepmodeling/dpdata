# Efficient LAMMPS Trajectory Frame Reading

This document describes the efficient trajectory frame reading functionality implemented for LAMMPS dump files in dpdata, addressing issue #367.

## Overview

The traditional approach to reading MD trajectories loads all frames into memory and then filters them. This can be inefficient when you only need specific frames from large trajectory files. The new implementation allows you to specify exactly which frames to read, skipping unwanted frames entirely.

## Key Features

### 1. Selective Frame Reading

Instead of loading entire trajectories, you can now specify exactly which frames to read:

```python
import dpdata

# Load only frames 23, 56, and 78 from a trajectory
system = dpdata.System(
    'trajectory.dump',
    fmt='lammps/dump',
    type_map=['O', 'H'],
    f_idx=[23, 56, 78]
)
```

### 2. Multi-Trajectory Pattern

The implementation supports the frames_dict pattern requested in the issue:

```python
import dpdata.lammps.dump as dump

frames_dict = {
    'trajectory1.dump': [23, 56, 78],
    'trajectory2.dump': [22],
    'trajectory3.dump': [10, 20, 30, 40]
}

# Load specified frames from multiple trajectories
data = dump.load_frames_from_trajectories(frames_dict, type_map=['O', 'H'])
```

### 3. Efficient Block Reading

The implementation uses block-based reading with `itertools.zip_longest` to skip frames efficiently:

- Determines frame structure (lines per frame) upfront
- Reads only requested frame blocks
- Skips unwanted frames without processing them

## API Reference

### Enhanced System Constructor

```python
dpdata.System(
    file_name,
    fmt='lammps/dump', 
    f_idx=None,  # NEW: List of frame indices to load
    **kwargs
)
```

**Parameters:**
- `f_idx` (list[int], optional): Specific frame indices to load (0-based). If provided, `begin` and `step` parameters are ignored.

### New Functions

#### `dpdata.lammps.dump.read_frames(fname, f_idx)`

Efficiently read specific frames from a LAMMPS dump file.

**Parameters:**
- `fname`: The dump file path
- `f_idx`: List of frame indices to read (0-based)

**Returns:**
- List of lines for the requested frames

#### `dpdata.lammps.dump.load_frames_from_trajectories(frames_dict, **kwargs)`

Load frames from multiple trajectory files using the frames_dict pattern.

**Parameters:**
- `frames_dict`: Dictionary mapping file paths to lists of frame indices
- `**kwargs`: Additional arguments passed to `system_data` (e.g., `type_map`, `unwrap`)

**Returns:**
- Combined system data dictionary

#### `dpdata.lammps.dump.get_frame_nlines(fname)`

Determine the number of lines per frame in a LAMMPS dump file.

**Parameters:**
- `fname`: The dump file path

**Returns:**
- Number of lines per frame (int)

## Performance Benefits

The efficient frame reading provides several advantages:

1. **Memory Efficiency**: Only loads requested frames into memory
2. **I/O Efficiency**: Skips unwanted frames during file reading
3. **Processing Efficiency**: No need to process and then discard unwanted frames

For large trajectory files with many frames, this can provide significant speedups when you only need a small subset of frames.

## Backward Compatibility

The implementation maintains full backward compatibility:

- Existing code using `begin` and `step` parameters continues to work unchanged
- All existing tests pass without modification
- The new `f_idx` parameter is optional and defaults to `None`

## Examples

### Basic Usage

```python
import dpdata

# Traditional approach (loads all frames)
system_all = dpdata.System('traj.dump', fmt='lammps/dump', type_map=['O', 'H'])

# Efficient approach (loads only specific frames)
system_subset = dpdata.System(
    'traj.dump', 
    fmt='lammps/dump', 
    type_map=['O', 'H'],
    f_idx=[10, 50, 100]
)
```

### Multi-Trajectory Loading

```python
import dpdata.lammps.dump as dump

# Define which frames to load from each trajectory
frames_dict = {
    'run1/traj.dump': [100, 200, 300],
    'run2/traj.dump': [50, 150, 250],
    'run3/traj.dump': [75, 175]
}

# Load all specified frames
data = dump.load_frames_from_trajectories(frames_dict, type_map=['C', 'H', 'O'])

# Convert to dpdata System if needed
system = dpdata.System(data=data)
```

### Performance Comparison

```python
import time
import dpdata

# Time traditional approach
start = time.time()
system = dpdata.System('large_traj.dump', fmt='lammps/dump', type_map=['O', 'H'])
filtered = system.sub_system([100, 500, 1000])
traditional_time = time.time() - start

# Time efficient approach
start = time.time()
system = dpdata.System(
    'large_traj.dump', 
    fmt='lammps/dump', 
    type_map=['O', 'H'],
    f_idx=[100, 500, 1000]
)
efficient_time = time.time() - start

print(f"Speedup: {traditional_time / efficient_time:.1f}x")
```

## Implementation Details

### Frame Structure Detection

The implementation first reads the file to determine the frame structure:

1. Finds the first "ITEM: TIMESTEP" line
2. Counts lines until the next "ITEM: TIMESTEP" 
3. Uses this count as the number of lines per frame

### Block-Based Reading

For selective frame reading:

1. Sorts requested frame indices for sequential access
2. Uses file position to skip to frame boundaries
3. Reads frame blocks only for requested indices
4. Combines results while preserving order

### Error Handling

The implementation handles various edge cases gracefully:

- Empty frame index lists return empty results
- Out-of-range indices are skipped silently
- Duplicate indices are automatically deduplicated
- Negative indices are ignored

This ensures robust operation even with invalid input.