#!/usr/bin/env python3
"""
Demonstration of efficient LAMMPS trajectory frame reading in dpdata.

This script shows how to use the new efficient frame reading functionality
that was implemented to address issue #367.
"""

import dpdata
import dpdata.lammps.dump as dump
import time


def demo_basic_usage():
    """Demonstrate basic usage of the new f_idx parameter."""
    print("=== Basic Usage Demo ===")
    
    # Traditional approach: load all frames
    print("1. Traditional approach - load all frames:")
    system_all = dpdata.System('tests/poscars/conf.5.dump', fmt='lammps/dump', type_map=['O', 'H'])
    print(f"   Loaded {len(system_all.data['coords'])} frames")
    
    # New efficient approach: load specific frames
    print("2. Efficient approach - load only frames [1, 3]:")
    system_selective = dpdata.System(
        'tests/poscars/conf.5.dump', 
        fmt='lammps/dump', 
        type_map=['O', 'H'], 
        f_idx=[1, 3]
    )
    print(f"   Loaded {len(system_selective.data['coords'])} frames")
    
    # Verify results are equivalent
    system_filtered = system_all.sub_system([1, 3])
    import numpy as np
    np.testing.assert_array_almost_equal(
        system_selective.data['coords'], 
        system_filtered.data['coords']
    )
    print("   ✓ Results match traditional filtering approach")


def demo_frames_dict_pattern():
    """Demonstrate the frames_dict pattern from the issue."""
    print("\n=== Frames Dict Pattern Demo ===")
    
    # This is the pattern requested in issue #367
    frames_dict = {
        'tests/poscars/conf.dump': [0, 1],        # Trajectory0: frames 0 and 1  
        'tests/poscars/conf.5.dump': [2, 4],      # Trajectory1: frames 2 and 4
    }
    
    print("Loading frames using the frames_dict pattern:")
    for traj, f_idx in frames_dict.items():
        print(f"  {traj}: frames {f_idx}")
    
    # Load using the new efficient function
    data = dump.load_frames_from_trajectories(frames_dict, type_map=['O', 'H'])
    
    print(f"Loaded {len(data['coords'])} frames total from {len(frames_dict)} trajectories")
    print("✓ Successfully combined frames from multiple trajectories")


def demo_performance_comparison():
    """Compare performance of different approaches."""
    print("\n=== Performance Comparison Demo ===")
    
    dump_file = 'tests/poscars/conf.5.dump'
    
    # Time the traditional approach
    start_time = time.time()
    system_all = dpdata.System(dump_file, fmt='lammps/dump', type_map=['O', 'H'])
    system_filtered = system_all.sub_system([1, 3])
    traditional_time = time.time() - start_time
    
    # Time the new efficient approach
    start_time = time.time()
    system_efficient = dpdata.System(
        dump_file, fmt='lammps/dump', type_map=['O', 'H'], f_idx=[1, 3]
    )
    efficient_time = time.time() - start_time
    
    print(f"Traditional (load all + filter): {traditional_time:.4f}s")
    print(f"Efficient (selective loading):   {efficient_time:.4f}s")
    
    if efficient_time < traditional_time:
        speedup = traditional_time / efficient_time
        print(f"✓ Speedup: {speedup:.1f}x faster")
    else:
        print("Note: For small files, the difference may not be noticeable")


def demo_api_usage():
    """Show various ways to use the new API."""
    print("\n=== API Usage Examples ===")
    
    # Method 1: Using dpdata.System with f_idx
    print("Method 1: dpdata.System with f_idx parameter")
    system = dpdata.System(
        'tests/poscars/conf.dump', 
        fmt='lammps/dump', 
        type_map=['O', 'H'], 
        f_idx=[1]
    )
    print(f"  Loaded {len(system.data['coords'])} frame(s)")
    
    # Method 2: Using the low-level read_frames function
    print("Method 2: Low-level read_frames function")
    lines = dump.read_frames('tests/poscars/conf.dump', [0, 1])
    data = dump.system_data(lines, type_map=['O', 'H'])
    print(f"  Loaded {len(data['coords'])} frame(s)")
    
    # Method 3: Using load_frames_from_trajectories for multiple files
    print("Method 3: load_frames_from_trajectories for multiple files")
    frames_dict = {'tests/poscars/conf.dump': [1]}
    data = dump.load_frames_from_trajectories(frames_dict, type_map=['O', 'H'])
    print(f"  Loaded {len(data['coords'])} frame(s)")


if __name__ == "__main__":
    print("LAMMPS Trajectory Efficient Frame Reading Demo")
    print("=" * 50)
    
    demo_basic_usage()
    demo_frames_dict_pattern() 
    demo_performance_comparison()
    demo_api_usage()
    
    print("\n" + "=" * 50)
    print("Demo completed! The new efficient frame reading functionality")
    print("allows you to load only the trajectory frames you need,")
    print("potentially saving significant time and memory for large files.")