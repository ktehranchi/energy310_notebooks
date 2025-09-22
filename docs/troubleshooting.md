# Troubleshooting Guide

## Common Issues and Solutions

### PyPSA Issues

#### Issue: "No solver found" error
**Symptoms**: Error when running `network.lopf()` or similar optimization commands.

**Solution**:
```bash
# Install CBC solver (open source)
conda install -c conda-forge coincbc

# Or install Gurobi (requires license)
conda install -c gurobi gurobi

# Test solver availability
python -c "import pypsa; print(pypsa.Network().solvers)"
```

#### Issue: Memory errors with large networks
**Symptoms**: Python crashes or "MemoryError" when solving large models.

**Solutions**:
1. Reduce temporal resolution:
   ```python
   # Instead of hourly, use 4-hourly snapshots
   network.set_snapshots(pd.date_range('2019', periods=2190, freq='4H'))
   ```

2. Reduce spatial detail:
   ```python
   # Aggregate buses or use clustered networks
   ```

3. Use linear approximation:
   ```python
   # Use LOPF instead of full AC power flow
   network.lopf()  # instead of network.pf()
   ```

#### Issue: Convergence problems
**Symptoms**: Solver runs but doesn't find a solution.

**Solutions**:
1. Check for infeasible constraints:
   ```python
   # Relax capacity constraints
   network.generators.loc[:, 'p_nom_extendable'] = True
   ```

2. Check data consistency:
   ```python
   # Ensure load equals generation potential
   total_load = network.loads_t.p_set.sum().sum()
   total_capacity = network.generators.p_nom.sum()
   print(f"Load: {total_load:.0f} MWh, Capacity: {total_capacity:.0f} MW")
   ```

### Git Issues

#### Issue: "Permission denied" when pushing
**Symptoms**: Cannot push to repository.

**Solutions**:
1. Check authentication:
   ```bash
   git config --global user.name "Your Name"
   git config --global user.email "your.email@example.com"
   ```

2. Use personal access token (not password) for GitHub

3. Check repository permissions

#### Issue: Merge conflicts
**Symptoms**: Git shows conflicts when merging or pulling.

**Solution**:
```bash
# Pull latest changes
git pull origin main

# Resolve conflicts in your editor
# Look for <<<<<<< ======= >>>>>>> markers

# Add resolved files
git add conflicted_file.py

# Complete the merge
git commit -m "Resolve merge conflicts"
```

### Snakemake Issues

#### Issue: "Rule not found" error
**Symptoms**: Snakemake cannot find specified rules.

**Solutions**:
1. Check file paths in Snakefile
2. Ensure you're in the correct directory
3. Verify rule names match exactly

#### Issue: "Missing input files" error
**Symptoms**: Snakemake cannot find input files.

**Solutions**:
1. Check file paths are correct
2. Run dry-run to see expected files:
   ```bash
   snakemake --dry-run
   ```
3. Create missing input files or adjust dependencies

#### Issue: Scripts fail within Snakemake
**Symptoms**: Individual scripts work but fail in Snakemake.

**Solutions**:
1. Check that scripts use `snakemake` object correctly:
   ```python
   # In script, access inputs and outputs via snakemake object
   input_file = snakemake.input[0]
   output_file = snakemake.output[0]
   ```

2. Ensure working directory is correct
3. Add error handling to scripts

### Jupyter Issues

#### Issue: Kernel dies or crashes
**Symptoms**: Jupyter notebook kernel restarts unexpectedly.

**Solutions**:
1. Reduce memory usage:
   ```python
   # Clear variables periodically
   del large_dataframe
   import gc; gc.collect()
   ```

2. Restart and run cells incrementally

3. Check for infinite loops or recursive calls

#### Issue: Plots not displaying
**Symptoms**: Matplotlib plots don't show in notebooks.

**Solutions**:
```python
# Ensure matplotlib backend is set
%matplotlib inline

# Or for interactive plots
%matplotlib widget
```

### Environment Issues

#### Issue: Package conflicts
**Symptoms**: Cannot install packages due to dependency conflicts.

**Solutions**:
1. Create fresh environment:
   ```bash
   conda deactivate
   conda env remove -n energy310
   conda env create -f environment.yml
   ```

2. Use mamba for faster solving:
   ```bash
   conda install mamba
   mamba env create -f environment.yml
   ```

#### Issue: Environment not found
**Symptoms**: Cannot activate conda environment.

**Solution**:
```bash
# List available environments
conda env list

# Recreate if missing
conda env create -f environment.yml
```

## Performance Optimization

### Large Model Handling

1. **Use efficient data types**:
   ```python
   # Use appropriate dtypes
   df['year'] = df['year'].astype('int16')  # instead of int64
   df['category'] = df['category'].astype('category')
   ```

2. **Lazy loading**:
   ```python
   # Load data only when needed
   def load_data(year):
       return pd.read_csv(f'data_{year}.csv')
   ```

3. **Parallel processing**:
   ```bash
   # Use multiple cores with Snakemake
   snakemake --cores 4
   ```

### Memory Management

1. **Monitor memory usage**:
   ```python
   import psutil
   print(f"Memory usage: {psutil.virtual_memory().percent}%")
   ```

2. **Clear unused variables**:
   ```python
   # In Jupyter notebooks
   %reset_selective -f "^temp_"  # Remove variables starting with temp_
   ```

## Debugging Tips

### PyPSA Debugging

1. **Check network consistency**:
   ```python
   network.consistency_check()
   ```

2. **Inspect components**:
   ```python
   print(network.buses)
   print(network.generators)
   print(network.loads)
   ```

3. **Check solver status**:
   ```python
   status = network.lopf()
   print(f"Solver status: {status}")
   ```

### General Python Debugging

1. **Use debugger**:
   ```python
   import pdb; pdb.set_trace()  # Set breakpoint
   ```

2. **Add logging**:
   ```python
   import logging
   logging.basicConfig(level=logging.INFO)
   logger = logging.getLogger(__name__)
   logger.info("Processing data...")
   ```

3. **Verbose output**:
   ```python
   # In scripts, add print statements
   print(f"Processing {len(data)} records...")
   ```

## Getting Additional Help

### Resources

- **PyPSA Documentation**: [https://pypsa.readthedocs.io/](https://pypsa.readthedocs.io/)
- **PyPSA Discussions**: [https://github.com/PyPSA/PyPSA/discussions](https://github.com/PyPSA/PyPSA/discussions)
- **Snakemake Documentation**: [https://snakemake.readthedocs.io/](https://snakemake.readthedocs.io/)
- **Stack Overflow**: Search for specific error messages

### Creating Bug Reports

When reporting issues, include:

1. **System information**:
   ```python
   import sys
   import pypsa
   print(f"Python: {sys.version}")
   print(f"PyPSA: {pypsa.__version__}")
   ```

2. **Minimal reproducible example**
3. **Full error message**
4. **Steps to reproduce the issue**

### Community Support

- GitHub Issues for bug reports
- GitHub Discussions for questions
- Course forum for class-specific questions

Remember: Most issues have been encountered before, so searching existing issues and documentation often provides quick solutions!