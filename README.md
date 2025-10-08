# Energy 310: Introduction to Energy Systems Modeling
### Stanford University, Energy Science and Engineering

Teaching Team:
 - Professor Inês Azevedo
 - Kamran Tehranchi, PhD Candidate


## 📚 Course Modules

1. **[What is PyPSA?](notebooks/01_what_is_pypsa.ipynb)** - Introduction to power system analysis with Python
2. **[Intro to PyPSA-USA]**
3. **[Introduction to snakemake]**
4 .... more on the way


## Installation

1. **Install uv (if not already installed):**
   ```bash
   # On macOS and Linux
   curl -LsSf https://astral.sh/uv/install.sh | sh
   ```

2. **Clone the repository:**
   ```bash
   git clone https://github.com/ktehranchi/energy310_notebooks.git
   cd energy310_notebooks
   ```

3. **Set up the environment using uv (recommended):**
   ```bash
   # UV automatically creates a virtual environment and installs dependencies
   uv sync

   # But because wwant to install extra dependancys specified in the pyproject use this command
   uv sync --extra geographic --extra advanced
   ```

4. **Open the Project in Your IDE or Jupyter Lab:**
   - **In Your IDE:** Open the project folder `energy310_notebooks` in your preferred IDE to start exploring the code.
   - **In Jupyter Lab:**
     ```bash
     # Using uv (runs in the project's virtual environment)
     uv run python -m ipykernel install --user --name energy310 --display-name "Python (energy310)"
     uv run jupyter lab
     ```

5. **Start with the first notebook:**
   Open `notebooks/01_what_is_pypsa.ipynb`

### What are UV Environments

UV creates **project-local virtual environments** that are automatically managed:

- **Automatic Environment Creation**: UV creates a `.venv` directory in your project folder
- **Dependency Management**: Dependencies are tracked in `pyproject.toml` and locked in `uv.lock`
- **No Manual Activation**: Use `uv run <command>` to execute commands in the environment
- **Fast Installation**: UV resolves and installs packages much faster than conda or pip
- **Reproducible Builds**: The `uv.lock` file ensures identical environments across machines

**Key UV commands:**
- `uv sync` - Install/update all dependencies
- `uv add <package>` - Add a new dependency
- `uv remove <package>` - Remove a dependency
- `uv run <command>` - Run a command in the project environment
- `uv shell` - Activate the environment in your current shell


## 🛠️ Key Technologies

- **[PyPSA](https://pypsa.org/)**: Open-source power system analysis framework
- **[Snakemake](https://snakemake.readthedocs.io/)**: Workflow management system
- **[Git](https://git-scm.com/)**: Version control system
- **Python Scientific Stack**: NumPy, Pandas, Matplotlib, Seaborn

## 📚 Additional Resources

### PyPSA Resources
- [PyPSA Documentation](https://pypsa.readthedocs.io/)
- [PyPSA-USA](https://github.com/PyPSA/pypsa-usa)
- [PyPSA Examples](https://github.com/PyPSA/pypsa-examples)
- [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur): Continental-scale model

### Git and Open Source
- [Git Documentation](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)

### Workflow Management
- [Snakemake Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)


## 📄 License

This course material is licensed under the MIT License. See [LICENSE](LICENSE) file for details.