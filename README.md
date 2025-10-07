# Energy 310: Electricity System Simulation and Optimization

Course materials for graduate-level energy systems analysis using PyPSA and modern software development practices.

## 📚 Course Overview

This repository contains three comprehensive modules designed to teach the fundamentals of electricity system modeling, optimization, and reproducible research workflows:

1. **[What is PyPSA?](notebooks/01_what_is_pypsa.ipynb)** - Introduction to power system analysis with Python
2. **[Git & Open Source Development](tutorials/02_git_and_opensource.md)** - Version control and collaborative research
3. **[Research Pipelines with Snakemake](notebooks/03_snakemake_pipelines.ipynb)** - Workflow management for reproducible analysis

## 🚀 Quick Start

### Installation

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

4. **Launch Jupyter Lab:**
   ```bash
   # Using uv (runs in the project's virtual environment)
   uv run jupyter lab
   ```

5. **Start with the first notebook:**
   Open `notebooks/01_what_is_pypsa.ipynb` and begin learning!

### Understanding UV Environments

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

## 📖 Module Descriptions

### Module 1: What is PyPSA?
**File:** `notebooks/01_what_is_pypsa.ipynb`

Learn the fundamentals of PyPSA (Python for Power System Analysis):
- Core components of power system models
- Building your first 3-bus network
- Running power flow analysis and optimization
- Time series analysis with renewable energy
- Visualization and results interpretation

**Learning Outcomes:**
- Understand power system modeling concepts
- Build and solve optimization problems
- Analyze renewable energy integration
- Create publication-quality visualizations

### Module 2: Git & Open Source Development
**File:** `tutorials/02_git_and_opensource.md`

Master version control and collaborative research:
- Open source software principles
- Git fundamentals and workflows
- Repository structure and best practices
- Contributing to open source projects
- Reproducible research practices

**Learning Outcomes:**
- Use Git for version control
- Collaborate on open source projects
- Structure research repositories
- Understand software licensing
- Implement reproducible workflows

### Module 3: Research Pipelines with Snakemake
**File:** `notebooks/03_snakemake_pipelines.ipynb`

Build scalable, reproducible analysis workflows:
- Workflow management principles
- Snakemake syntax and concepts
- Energy system analysis pipelines
- Parallel computation and scaling
- Integration with PyPSA

**Learning Outcomes:**
- Design modular analysis workflows
- Parallelize computations
- Manage complex dependencies
- Scale from laptop to cluster
- Ensure reproducible results

## 🗂️ Repository Structure

```
energy310_notebooks/
├── README.md                           # This file
├── pyproject.toml                      # UV project configuration and dependencies
├── environment.yml                     # Legacy conda environment specification
├── requirements.txt                    # Python package requirements
├── .gitignore                         # Git ignore file
│
├── notebooks/                         # Jupyter notebooks
│   ├── 01_what_is_pypsa.ipynb        # PyPSA introduction
│   └── 03_snakemake_pipelines.ipynb  # Workflow management
│
├── tutorials/                         # Tutorial materials
│   └── 02_git_and_opensource.md      # Git and open source guide
│
├── workflows/                         # Snakemake workflows
│   ├── Snakefile                      # Main workflow
│   ├── config.yaml                    # Configuration file
│   └── scripts/                       # Analysis scripts
│       └── download_data.py           # Data preparation script
│
├── data/                              # Data directory
│   ├── examples/                      # Example datasets
│   ├── raw/                          # Raw input data
│   └── processed/                     # Processed data
│
├── docs/                              # Additional documentation
└── results/                           # Analysis outputs
    ├── plots/                         # Generated figures
    └── tables/                        # Summary tables
```

## 🛠️ Key Technologies

- **[PyPSA](https://pypsa.org/)**: Open-source power system analysis framework
- **[Snakemake](https://snakemake.readthedocs.io/)**: Workflow management system
- **[Jupyter](https://jupyter.org/)**: Interactive computing environment
- **[Git](https://git-scm.com/)**: Version control system
- **Python Scientific Stack**: NumPy, Pandas, Matplotlib, Seaborn

## 📋 Learning Path

### Beginner Path
1. Start with PyPSA notebook to understand energy modeling
2. Learn Git basics for version control
3. Build simple Snakemake workflows

### Advanced Path
1. Complete all modules in order
2. Implement your own research project
3. Contribute to open source energy tools
4. Deploy workflows on high-performance computing systems

## 🎯 Assignments and Exercises

Each module includes practical exercises:

- **Module 1**: Build progressively complex PyPSA networks
- **Module 2**: Set up a research repository with proper Git workflow
- **Module 3**: Create a complete analysis pipeline for your research

## 📚 Additional Resources

### PyPSA Resources
- [PyPSA Documentation](https://pypsa.readthedocs.io/)
- [PyPSA Examples](https://github.com/PyPSA/pypsa-examples)
- [PyPSA-Eur](https://github.com/PyPSA/pypsa-eur): Continental-scale model

### Git and Open Source
- [Git Documentation](https://git-scm.com/doc)
- [GitHub Guides](https://guides.github.com/)
- [The Turing Way](https://the-turing-way.netlify.app/): Reproducible research guide

### Workflow Management
- [Snakemake Tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html)
- [Best Practices for Scientific Computing](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.1001745)

## 🤝 Contributing

We welcome contributions! Please see our contributing guidelines:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

For bug reports and feature requests, please use GitHub Issues.

## 📄 License

This course material is licensed under the MIT License. See [LICENSE](LICENSE) file for details.

## 👥 Authors and Acknowledgments

- **Course Development**: Energy 310 Teaching Team
- **PyPSA Community**: For the excellent power system analysis framework
- **Open Source Contributors**: For tools and inspiration

## 🆘 Getting Help

- **Course Questions**: Use GitHub Issues or discussions
- **Technical Issues**: Check the troubleshooting section in each module
- **PyPSA Help**: [PyPSA Discussion Forum](https://github.com/PyPSA/PyPSA/discussions)
- **Git Help**: [Git Community](https://git-scm.com/community)

## 🔄 Course Updates

This repository is actively maintained. Check the [Releases](https://github.com/ktehranchi/energy310_notebooks/releases) page for updates and new content.

---

**Happy Learning! 🎓⚡**

*Empowering the next generation of energy system researchers with modern computational tools and practices.*