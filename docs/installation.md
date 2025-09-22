# Installation Guide

## System Requirements

- **Operating System**: Windows 10+, macOS 10.15+, or Linux
- **Python**: 3.9 or higher
- **Memory**: Minimum 8GB RAM (16GB recommended for large models)
- **Storage**: At least 5GB free space

## Installation Methods

### Method 1: Conda (Recommended)

Conda automatically handles complex dependencies and is the most reliable installation method.

```bash
# Clone the repository
git clone https://github.com/ktehranchi/energy310_notebooks.git
cd energy310_notebooks

# Create and activate environment
conda env create -f environment.yml
conda activate energy310

# Verify installation
python -c "import pypsa; print(f'PyPSA version: {pypsa.__version__}')"
```

### Method 2: pip

If you prefer pip or don't have conda:

```bash
# Clone the repository
git clone https://github.com/ktehranchi/energy310_notebooks.git
cd energy310_notebooks

# Create virtual environment (optional but recommended)
python -m venv energy310_env
source energy310_env/bin/activate  # On Windows: energy310_env\Scripts\activate

# Install requirements
pip install -r requirements.txt

# Verify installation
python -c "import pypsa; print(f'PyPSA version: {pypsa.__version__}')"
```

## Optional Dependencies

### Optimization Solvers

For better performance, install commercial solvers:

- **Gurobi** (free academic license): [Gurobi Installation](https://www.gurobi.com/academia/academic-program-and-licenses/)
- **CPLEX** (free academic license): [CPLEX Installation](https://www.ibm.com/academic/home)

If you don't have commercial solvers, the open-source CBC solver is included.

### Geographic Visualization

For map plotting capabilities:

```bash
conda install -c conda-forge cartopy geopandas
```

## Troubleshooting

### Common Issues

#### 1. Solver Not Found
```
Error: No solver found
```

**Solution**: Install CBC solver:
```bash
conda install -c conda-forge coincbc
```

#### 2. Import Error for PyPSA
```
ModuleNotFoundError: No module named 'pypsa'
```

**Solution**: Ensure you've activated the correct environment:
```bash
conda activate energy310
```

#### 3. Memory Issues
If you encounter memory errors with large models:
- Increase virtual memory/swap space
- Use a machine with more RAM
- Reduce model complexity

#### 4. Jupyter Lab Not Starting
```bash
# If jupyter lab command not found
conda install -c conda-forge jupyterlab

# Or
pip install jupyterlab
```

### Platform-Specific Issues

#### Windows
- Use Anaconda Prompt instead of Command Prompt
- Some packages might require Microsoft Visual C++ Redistributable

#### macOS
- You might need to install Xcode Command Line Tools:
  ```bash
  xcode-select --install
  ```

#### Linux
- Ensure build tools are installed:
  ```bash
  sudo apt-get install build-essential  # Ubuntu/Debian
  sudo yum groupinstall "Development Tools"  # CentOS/RHEL
  ```

## Testing Your Installation

Run this test script to verify everything works:

```python
# test_installation.py
import sys
print(f"Python version: {sys.version}")

try:
    import pypsa
    print(f"✓ PyPSA {pypsa.__version__} installed successfully")
except ImportError:
    print("✗ PyPSA not found")

try:
    import pandas as pd
    print(f"✓ Pandas {pd.__version__} installed successfully")
except ImportError:
    print("✗ Pandas not found")

try:
    import matplotlib
    print(f"✓ Matplotlib {matplotlib.__version__} installed successfully")
except ImportError:
    print("✗ Matplotlib not found")

try:
    import snakemake
    print(f"✓ Snakemake {snakemake.__version__} installed successfully")
except ImportError:
    print("✗ Snakemake not found")

print("\nTesting PyPSA functionality...")
try:
    network = pypsa.Network()
    network.add("Bus", "test_bus")
    print("✓ PyPSA basic functionality works")
except Exception as e:
    print(f"✗ PyPSA test failed: {e}")

print("\nInstallation test complete!")
```

## Getting Help

If you encounter issues not covered here:

1. Check the [GitHub Issues](https://github.com/ktehranchi/energy310_notebooks/issues)
2. Search the [PyPSA Documentation](https://pypsa.readthedocs.io/)
3. Ask on the [PyPSA Discussion Forum](https://github.com/PyPSA/PyPSA/discussions)
4. Create a new issue with your system details and error messages