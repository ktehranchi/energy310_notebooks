# Introduction to Open-Source Software Development and Git

## Learning Objectives

By the end of this tutorial, you will:
1. Understand the principles of open-source software development
2. Learn the basics of version control with Git
3. Know how to collaborate on open-source projects
4. Understand the importance of reproducible research
5. Learn best practices for scientific computing with version control

## 1. What is Open-Source Software?

Open-source software (OSS) is software whose source code is freely available for anyone to view, modify, and distribute. This model has become the foundation of modern scientific computing and research.

### Key Principles of Open Source:

- **Transparency**: Source code is open for inspection
- **Collaboration**: Multiple developers can contribute
- **Community**: Users become contributors and maintainers
- **Reproducibility**: Others can verify and build upon research
- **Innovation**: Rapid development through shared knowledge

### Examples in Energy Research:

- **PyPSA**: Power system analysis and optimization
- **HOMER**: Hybrid renewable energy system design
- **EnergyPLAN**: Energy system analysis
- **Calliope**: Multi-scale energy system modeling
- **OSeMOSYS**: Long-term energy system optimization

## 2. Why Open Source Matters in Academic Research

### Reproducibility Crisis
Modern research faces challenges in reproducibility. Open-source tools and transparent methodologies help address this by:
- Making methods publicly available
- Allowing others to verify results
- Enabling building upon previous work
- Reducing "black box" approaches

### Benefits for Researchers:
- **Faster Development**: Build on existing tools
- **Quality Improvement**: Community review and testing
- **Career Benefits**: Public contributions visible to employers
- **Learning Opportunities**: Read and learn from expert code
- **Citation Impact**: Open tools get cited more often

## 3. Introduction to Git

Git is a distributed version control system that tracks changes in files over time. It's essential for:
- Managing code versions
- Collaborating with others
- Maintaining project history
- Recovering from mistakes

### Key Git Concepts:

#### Repository (Repo)
A project folder tracked by Git, containing all files and their history.

#### Commit
A snapshot of your project at a specific point in time, with a descriptive message.

#### Branch
A parallel version of your repository where you can make changes without affecting the main code.

#### Merge
Combining changes from one branch into another.

#### Remote
A version of your repository hosted on a server (like GitHub, GitLab).

## 4. Git Workflow for Research Projects

### Basic Git Commands

```bash
# Initialize a new repository
git init

# Clone an existing repository
git clone https://github.com/username/repository.git

# Check status of your files
git status

# Add files to staging area
git add filename.py
git add .  # Add all files

# Commit changes
git commit -m "Descriptive message about changes"

# Push changes to remote repository
git push origin main

# Pull latest changes from remote
git pull origin main

# View commit history
git log --oneline
```

### Research-Specific Workflow

1. **Create a new branch for each experiment/feature**
   ```bash
   git checkout -b experiment-wind-storage
   ```

2. **Make incremental commits with clear messages**
   ```bash
   git add analysis.py
   git commit -m "Add wind-storage optimization model"
   ```

3. **Push branch and create pull request for review**
   ```bash
   git push origin experiment-wind-storage
   ```

4. **Merge successful experiments into main branch**

## 5. Best Practices for Research with Git

### Commit Messages
Write clear, descriptive commit messages:

**Good examples:**
- `Add renewable energy potential analysis for California`
- `Fix bug in power flow calculation for DC lines`
- `Update documentation for storage optimization`

**Poor examples:**
- `Fixed stuff`
- `Update`
- `Changes`

### Repository Structure
Organize your research repository logically:

```
research-project/
├── README.md
├── environment.yml
├── requirements.txt
├── data/
│   ├── raw/
│   ├── processed/
│   └── external/
├── notebooks/
│   ├── 01_data_exploration.ipynb
│   ├── 02_model_development.ipynb
│   └── 03_results_analysis.ipynb
├── src/
│   ├── __init__.py
│   ├── data_processing.py
│   ├── models.py
│   └── visualization.py
├── tests/
│   ├── test_data_processing.py
│   └── test_models.py
├── docs/
│   └── methodology.md
└── results/
    ├── figures/
    └── tables/
```

### What to Track and What to Ignore

**Track:**
- Source code (.py, .R, .m files)
- Documentation (.md, .rst files)
- Configuration files
- Small data files
- Environment specifications

**Don't Track (use .gitignore):**
- Large data files (use Git LFS or external storage)
- Temporary files
- Log files
- Personal IDE settings
- Compiled files

## 6. GitHub for Research Collaboration

### Key Features:

#### Issues
Track bugs, feature requests, and discussions
- Use labels to categorize (bug, enhancement, question)
- Assign to team members
- Reference in commits with #issue-number

#### Pull Requests
Propose and review changes before merging
- Describe what changes and why
- Request reviews from collaborators
- Use automated testing

#### Releases
Tag stable versions of your research
- Use semantic versioning (v1.0.0, v1.1.0, etc.)
- Include release notes
- Attach compiled results or data

#### Documentation
Use README files, wikis, and GitHub Pages
- Installation instructions
- Usage examples
- Citation information
- Contributing guidelines

## 7. Open Source Licenses for Research

### Common Licenses:

#### MIT License
- Very permissive
- Allows commercial use
- Minimal restrictions

#### Apache 2.0
- Permissive with patent protection
- Good for larger projects
- Widely used in industry

#### GPL (GNU General Public License)
- Copyleft license
- Derivative works must also be open source
- Strong community protection

#### Creative Commons (for data/documents)
- Various levels of openness
- Common for research data and publications

### Choosing a License:
- Consider your institution's policies
- Think about how you want others to use your work
- Check compatibility with dependencies
- When in doubt, consult legal counsel

## 8. Contributing to Open Source Projects

### How to Get Started:

1. **Find Projects**: Look for projects you use in your research
2. **Start Small**: Fix typos, improve documentation
3. **Read Contributing Guidelines**: Each project has different processes
4. **Follow Code Style**: Maintain consistency with existing code
5. **Test Your Changes**: Ensure you don't break existing functionality

### Example Contribution Workflow:

```bash
# 1. Fork the repository on GitHub
# 2. Clone your fork
git clone https://github.com/yourusername/pypsa.git

# 3. Create a new branch
git checkout -b fix-documentation-typo

# 4. Make your changes
# Edit files...

# 5. Test your changes
python -m pytest tests/

# 6. Commit and push
git add .
git commit -m "Fix typo in renewable energy documentation"
git push origin fix-documentation-typo

# 7. Create pull request on GitHub
```

## 9. Version Control for Data and Results

### Challenges with Research Data:
- Large file sizes
- Binary formats
- Frequent changes
- Multiple versions

### Solutions:

#### Git LFS (Large File Storage)
For large files that need version control:
```bash
git lfs install
git lfs track "*.csv"
git add .gitattributes
git add large_dataset.csv
git commit -m "Add large dataset with LFS"
```

#### External Storage with Versioning
- AWS S3 with versioning
- Google Cloud Storage
- Institutional data repositories
- Zenodo for permanent storage

#### Data Version Control (DVC)
Specialized tool for data science:
```bash
dvc add data/large_dataset.csv
git add data/large_dataset.csv.dv .gitignore
git commit -m "Add dataset with DVC"
```

## 10. Reproducible Research Practices

### Environment Management

#### UV Project Configuration (Recommended)
Modern Python package management with fast dependency resolution:

```toml
# pyproject.toml
[project]
name = "energy-research"
version = "0.1.0"
requires-python = ">=3.9"
dependencies = [
    "pypsa>=0.21.0",
    "pandas>=1.3.0,<2.0.0",
    "matplotlib>=3.4.0",
    "seaborn>=0.11.0",
]

[tool.uv]
dev-dependencies = [
    "pytest",
    "black",
    "jupyter",
]
```

#### Conda Environment Files (Legacy)
```yaml
name: energy-research
channels:
  - conda-forge
dependencies:
  - python=3.9
  - pypsa=0.21.0
  - pandas=1.3.0
  - matplotlib=3.4.0
  - pip:
    - custom-package==1.0.0
```

#### Requirements Files
```txt
pypsa==0.21.0
pandas>=1.3.0,<2.0.0
matplotlib>=3.4.0
seaborn>=0.11.0
```

#### Environment Setup Commands

```bash
# Using UV (recommended)
uv sync                    # Install/sync dependencies
uv add new-package        # Add a new dependency
uv run python script.py   # Run in environment
uv shell                  # Activate shell

# Using conda (legacy)
conda env create -f environment.yml
conda activate energy-research

# Using pip
pip install -r requirements.txt
```

### Computational Notebooks

#### Best Practices:
- Clear markdown explanations
- Reproducible execution order
- Version control notebooks (consider nbstripout)
- Include requirements and environment info

#### Jupyter Extensions:
- nbdime for better diffs
- nbstripout to clean notebooks before committing
- papermill for parameterized notebook execution

### Documentation

#### README Template:
```markdown
# Project Title

Brief description of your research project.

## Installation

### Requirements
- Python 3.9+
- See pyproject.toml for packages

### Setup using UV (recommended)
```bash
git clone https://github.com/username/project.git
cd project
uv sync  # Creates environment and installs dependencies
```

### Alternative setup methods
```bash
# Using conda (legacy)
conda env create -f environment.yml
conda activate energy-research

# Using pip
pip install -r requirements.txt
```

## Usage

```bash
# Run scripts with UV
uv run python analysis.py

# Start Jupyter
uv run jupyter lab

# Or activate environment first
uv shell
python analysis.py
```

## Usage

Basic usage examples...

## Data

Description of data sources and how to obtain them...

## Results

How to reproduce key results...

## Citation

How to cite this work...

## License

MIT License - see LICENSE file
```

## 11. Practical Exercise

### Setting up Your First Research Repository

1. **Create a new repository on GitHub**
   - Include README, .gitignore (Python), and license

2. **Clone it locally**
   ```bash
   git clone https://github.com/yourusername/energy-research-template.git
   cd energy-research-template
   ```

3. **Set up the directory structure**
   ```bash
   mkdir -p data/{raw,processed} notebooks src tests docs results
   touch src/__init__.py
   ```

4. **Set up project configuration**
   - Create a `pyproject.toml` file with your dependencies
   - Or copy the `pyproject.toml` from this course repository
   - Initialize the UV project: `uv sync`

5. **Write a proper README**
   - Describe your research question
   - Include installation instructions
   - Add usage examples

6. **Create your first notebook**
   - Add it to the notebooks/ directory
   - Include clear markdown documentation

7. **Commit and push everything**
   ```bash
   git add .
   git commit -m "Initial repository setup for energy research"
   git push origin main
   ```

## 12. Advanced Topics

### Branching Strategies

#### Feature Branch Workflow
- Create branches for each new feature/experiment
- Keep main branch stable
- Use pull requests for code review

#### Git Flow
- More complex branching model
- Separate branches for development, features, releases
- Good for larger projects with multiple collaborators

### Automation and CI/CD

#### GitHub Actions
Automate testing and deployment:

```yaml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    steps:
    - uses: actions/checkout@v2
    - name: Set up Python
      uses: actions/setup-python@v2
      with:
        python-version: 3.9
    - name: Install dependencies
      run: |
        pip install -r requirements.txt
    - name: Run tests
      run: |
        python -m pytest tests/
```

### Collaboration Tools

#### Code Review
- Use pull requests for all changes
- Provide constructive feedback
- Focus on readability and reproducibility

#### Project Management
- Use GitHub Issues for task tracking
- Create milestones for project phases
- Use project boards for workflow management

## 13. Summary and Best Practices

### Key Takeaways:

1. **Start Early**: Begin using version control from day one of your research
2. **Commit Often**: Small, frequent commits are better than large, infrequent ones
3. **Document Everything**: Good documentation is as important as good code
4. **Test Your Work**: Automated testing prevents errors and increases confidence
5. **Share Your Work**: Open source contributions benefit your career and the community
6. **Learn Continuously**: Version control and open source practices are constantly evolving

### Research-Specific Guidelines:

- **Track methodology changes**: Use branches for different approaches
- **Version your data**: Keep track of data processing steps
- **Archive releases**: Create releases for paper submissions and revisions
- **Cite software**: Properly cite the open-source tools you use
- **Give back**: Contribute improvements to tools you use

### Common Pitfalls to Avoid:

- Committing large binary files without LFS
- Writing vague commit messages
- Not backing up your work to remote repositories
- Ignoring merge conflicts
- Not reading project contribution guidelines

### Resources for Continued Learning:

- **Git Documentation**: [https://git-scm.com/doc](https://git-scm.com/doc)
- **GitHub Guides**: [https://guides.github.com/](https://guides.github.com/)
- **Open Source Guide**: [https://opensource.guide/](https://opensource.guide/)
- **The Turing Way**: Guide to reproducible research
- **Software Carpentry**: Workshops on research computing skills

## 14. Assignment Ideas

1. **Repository Setup**: Create and organize a research repository with proper structure
2. **Contribution Practice**: Find and fix a small issue in an open-source energy tool
3. **Collaboration Exercise**: Work in pairs to merge conflicting changes
4. **Documentation Challenge**: Write comprehensive documentation for a code snippet
5. **License Analysis**: Compare different licenses for a hypothetical research project

Remember: Good version control practices and open source collaboration skills are invaluable for modern research careers. Start building these habits early!