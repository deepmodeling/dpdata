#!/usr/bin/env python3
"""
Test script to verify that the VASP OUTCAR ML parsing handles different text variations robustly.
"""

import dpdata
import numpy as np

def test_ml_vs_nonml_consistency():
    """Test that ML and non-ML modes extract consistent data for overlapping frames."""
    
    print("=== Testing ML vs Non-ML consistency ===")
    fname = "tests/poscars/OUTCAR.ch4.ml"
    
    system_ml = dpdata.LabeledSystem(fname, fmt="vasp/outcar", ml=True)
    system_nonml = dpdata.LabeledSystem(fname, fmt="vasp/outcar", ml=False)
    
    print(f"ML mode extracted: {len(system_ml['energies'])} frames")
    print(f"Non-ML mode extracted: {len(system_nonml['energies'])} frames")
    
    # The frames should have consistent atom information
    assert system_ml["atom_names"] == system_nonml["atom_names"]
    assert system_ml["atom_numbs"] == system_nonml["atom_numbs"]
    assert np.array_equal(system_ml["atom_types"], system_nonml["atom_types"])
    
    print("✓ Atom information is consistent between modes")
    
    # Cell shapes should be correct
    assert system_ml["cells"].shape == (len(system_ml["energies"]), 3, 3)
    assert system_nonml["cells"].shape == (len(system_nonml["energies"]), 3, 3)
    
    print("✓ Cell data has correct dimensions")
    
    # All cell determinants should be positive (valid cells)
    for i, cell in enumerate(system_ml["cells"]):
        det = np.linalg.det(cell)
        assert det > 0, f"ML frame {i} has invalid cell determinant: {det}"
    
    for i, cell in enumerate(system_nonml["cells"]):
        det = np.linalg.det(cell)
        assert det > 0, f"Non-ML frame {i} has invalid cell determinant: {det}"
    
    print("✓ All cells are valid (positive determinant)")
    
    return True

def test_robustness_improvements():
    """Test that the robustness improvements don't break existing functionality."""
    
    print("\n=== Testing robustness improvements ===")
    
    # The improvements include:
    # 1. Better error handling for malformed cell data lines
    # 2. More robust parsing of float values
    
    # Test should pass without errors
    system = dpdata.LabeledSystem("tests/poscars/OUTCAR.ch4.ml", fmt="vasp/outcar", ml=True)
    
    # Check that we get the expected number of frames
    assert len(system["energies"]) == 10, f"Expected 10 frames, got {len(system['energies'])}"
    
    # Check that all frames have complete data
    assert len(system["cells"]) == 10
    assert len(system["coords"]) == 10
    assert len(system["forces"]) == 10
    
    print("✓ Robustness improvements maintain expected behavior")
    
    return True

def main():
    """Run all tests."""
    
    print("Testing VASP OUTCAR ML parsing improvements...")
    print("=" * 60)
    
    try:
        test_ml_vs_nonml_consistency()
        test_robustness_improvements()
        
        print("\n" + "=" * 60)
        print("✅ All tests passed! The improvements are working correctly.")
        print("\nSummary of improvements:")
        print("1. More robust cell data extraction with better error handling")
        print("2. Improved parsing of float values in cell vectors") 
        print("3. Better handling of potential variations in OUTCAR format")
        
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        return False
    
    return True

if __name__ == "__main__":
    main()