# Installation Guide

## System Requirements

- **Operating System**: Windows 10+, macOS 10.15+, or Linux
- **Python**: 3.9 or higher
- **Memory**: Minimum 8GB RAM (16GB recommended for large models)
- **Storage**: At least 5GB free space

## Installation Methods

### Method 1: UV (Recommended)

UV is an extremely fast Python package manager that automatically handles dependencies and virtual environments.

#### Installing UV

```bash
# On macOS and Linux
curl -LsSf https://astral.sh/uv/install.sh | sh

# On Windows (PowerShell)
powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# Or using pip (if you already have Python)
pip install uv
```

#### Setting up the project

```bash
# Clone the repository
git clone https://github.com/ktehranchi/energy310_notebooks.git
cd energy310_notebooks

# Create environment and install dependencies
uv sync

# Verify installation
uv run python -c "import pypsa; print(f'PyPSA version: {pypsa.__version__}')"
```

### Method 2: pip

If you prefer pip or don't want to use uv:

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

### Legacy Method: Conda

For compatibility, conda installation is still supported:

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

## Understanding UV Environments

UV provides a modern approach to Python environment management:

### How UV Environments Work

1. **Project-Local Environments**: UV creates a `.venv` directory within your project folder
2. **Automatic Management**: No need to manually activate/deactivate environments
3. **Dependency Tracking**: Dependencies are specified in `pyproject.toml` and locked in `uv.lock`
4. **Reproducible Builds**: The lock file ensures identical environments across different machines
5. **Fast Resolution**: UV's Rust-based resolver is significantly faster than pip or conda

### Key UV Commands

```bash
# Install/sync all dependencies
uv sync

# Add a new dependency
uv add numpy matplotlib

# Remove a dependency
uv remove old-package

# Run commands in the environment
uv run python script.py
uv run jupyter lab

# Activate shell in the environment (optional)
uv shell

# Update dependencies to latest compatible versions
uv sync --upgrade

# Export to requirements.txt format
uv export --format requirements-txt --output-file requirements.txt
```

### UV vs Conda/Pip

| Feature | UV | Conda | Pip |
|---------|----|----|-----|
| Speed | ⚡ Very Fast | 🐌 Slow | 🚶 Medium |
| Environment Management | ✅ Automatic | ⚠️ Manual | ⚠️ Manual |
| Dependency Resolution | ✅ Advanced | ✅ Good | ⚠️ Basic |
| Cross-platform | ✅ Yes | ✅ Yes | ✅ Yes |
| Binary packages | ⚠️ PyPI only | ✅ Multiple channels | ⚠️ PyPI only |

## Optional Dependencies

### Optimization Solvers

For better performance, install commercial solvers:

- **Gurobi** (free academic license): [Gurobi Installation](https://www.gurobi.com/academia/academic-program-and-licenses/)
- **CPLEX** (free academic license): [CPLEX Installation](https://www.ibm.com/academic/home)

If you don't have commercial solvers, the open-source CBC solver is included.

### Geographic Visualization

For map plotting capabilities:

```bash
# Using UV
uv add cartopy geopandas

# Using conda (if using legacy method)
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
# Using UV
uv add pulp  # CBC solver is included with PuLP

# Or install additional solvers (if available via pip)
uv add gurobipy  # Requires Gurobi license

# Using conda (legacy method)
conda install -c conda-forge coincbc
```

#### 2. Import Error for PyPSA
```
ModuleNotFoundError: No module named 'pypsa'
```

**Solution**: Ensure you're using the project environment:
```bash
# Using UV - run in project environment
uv run python -c "import pypsa"

# Or activate the environment
uv shell
python -c "import pypsa"

# Using conda (legacy method)
conda activate energy310
```

#### 3. Memory Issues
If you encounter memory errors with large models:
- Increase virtual memory/swap space
- Use a machine with more RAM
- Reduce model complexity

#### 4. Jupyter Lab Not Starting
```bash
# If jupyter lab command not found with UV
uv add jupyterlab
uv run jupyter lab

# Or using pip
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

Run this verification script to ensure everything works correctly:

```bash
# Using UV
uv run python verify_installation.py

# Using conda (legacy method)
python verify_installation.py
```

The script will test all core packages and PyPSA functionality. You should see ✓ marks for all core packages.

### Manual Testing

You can also test manually:

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