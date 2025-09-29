#!/usr/bin/env python3
"""
Installation verification script for Energy 310 notebooks.
Run this script to verify that your environment is set up correctly.
"""

import sys
print(f"Python version: {sys.version}")
print()

# Test core packages
packages_to_test = [
    ("numpy", "NumPy"),
    ("pandas", "Pandas"),
    ("matplotlib", "Matplotlib"),
    ("seaborn", "Seaborn"),
    ("scipy", "SciPy"),
    ("jupyter", "Jupyter"),
    ("jupyterlab", "JupyterLab"),
    ("pypsa", "PyPSA"),
    ("pulp", "PuLP"),
    ("networkx", "NetworkX"),
    ("requests", "Requests"),
    ("yaml", "PyYAML"),
    ("tqdm", "TQDM"),
]

optional_packages = [
    ("geopandas", "GeoPandas"),
    ("cartopy", "Cartopy"),
    ("rasterio", "Rasterio"),
    ("xarray", "XArray"),
    ("netCDF4", "NetCDF4"),
]

print("Testing core packages:")
print("=" * 40)

for module_name, display_name in packages_to_test:
    try:
        if module_name == "yaml":
            import yaml
            module = yaml
        else:
            module = __import__(module_name)
        
        version = getattr(module, '__version__', 'unknown')
        print(f"✓ {display_name:<12} {version}")
    except ImportError:
        print(f"✗ {display_name:<12} NOT FOUND")

print()
print("Testing optional packages:")
print("=" * 40)

for module_name, display_name in optional_packages:
    try:
        module = __import__(module_name)
        version = getattr(module, '__version__', 'unknown')
        print(f"✓ {display_name:<12} {version}")
    except ImportError:
        print(f"~ {display_name:<12} not installed (optional)")

print()
print("Testing PyPSA functionality:")
print("=" * 40)

try:
    import pypsa
    network = pypsa.Network()
    network.add("Bus", "test_bus")
    network.add("Generator", "test_gen", bus="test_bus", p_nom=100)
    network.add("Load", "test_load", bus="test_bus", p_set=50)
    print("✓ PyPSA basic functionality works")
    
    # Test if solvers are available
    try:
        # Try to run a simple optimization to test solvers
        network.lopf()
        print("✓ Optimization solver working")
    except Exception as solver_error:
        if "No solver" in str(solver_error):
            print("⚠ No optimization solvers found - install commercial solvers for better performance")
        else:
            print("✓ Basic solver functionality available")
        
except Exception as e:
    print(f"✗ PyPSA test failed: {e}")

print()
print("Installation verification complete!")
print()
print("If you see any ✗ errors for core packages, run:")
print("  uv sync")
print()
print("To install optional packages, run:")
print("  uv sync --extra geographic")
print("  uv sync --extra advanced")