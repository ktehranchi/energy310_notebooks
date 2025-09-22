#!/usr/bin/env python3
"""
Download and prepare input data for energy system analysis
"""

import pandas as pd
import numpy as np
import requests
import os
from pathlib import Path

def create_sample_demand_data():
    """Create sample demand data for demonstration"""
    # Create hourly demand data for a year
    hours = pd.date_range('2019-01-01', periods=8760, freq='H')
    
    # Create realistic demand profile
    # Base load + daily pattern + seasonal pattern + noise
    base_load = 5000  # MW
    
    # Daily pattern (peak at 18:00, minimum at 4:00)
    daily_pattern = 0.3 * np.sin(2 * np.pi * (hours.hour - 6) / 24)
    
    # Seasonal pattern (higher in summer and winter)
    seasonal_pattern = 0.2 * np.sin(2 * np.pi * (hours.dayofyear - 81) / 365)
    
    # Weekly pattern (lower on weekends)
    weekly_pattern = np.where(hours.weekday < 5, 0.1, -0.05)
    
    # Random noise
    noise = 0.05 * np.random.randn(len(hours))
    
    # Combine all patterns
    demand = base_load * (1 + daily_pattern + seasonal_pattern + weekly_pattern + noise)
    
    # Ensure positive demand
    demand = np.maximum(demand, base_load * 0.3)
    
    df = pd.DataFrame({
        'timestamp': hours,
        'demand_mw': demand,
        'region': 'central'
    })
    
    return df

def create_sample_renewable_data():
    """Create sample renewable generation profiles"""
    hours = pd.date_range('2019-01-01', periods=8760, freq='H')
    
    # Solar profile (sine wave during day, zero at night)
    solar_base = np.maximum(0, np.sin(2 * np.pi * (hours.hour - 6) / 24))
    
    # Add seasonal variation (higher in summer)
    solar_seasonal = 1 + 0.3 * np.sin(2 * np.pi * (hours.dayofyear - 81) / 365)
    
    # Add weather variability
    solar_weather = np.random.beta(2, 2, len(hours))
    
    solar_cf = solar_base * solar_seasonal * solar_weather * 0.25  # Scale to realistic capacity factors
    solar_cf = np.clip(solar_cf, 0, 1)
    
    # Wind profile (more variable, some seasonal pattern)
    wind_base = 0.3 + 0.2 * np.sin(2 * np.pi * hours.dayofyear / 365)  # Higher in winter
    wind_variability = np.random.weibull(2, len(hours)) * 0.5
    wind_cf = wind_base + wind_variability
    wind_cf = np.clip(wind_cf, 0, 1)
    
    df = pd.DataFrame({
        'timestamp': hours,
        'solar_cf': solar_cf,
        'wind_onshore_cf': wind_cf,
        'wind_offshore_cf': wind_cf * 1.3,  # Higher offshore capacity factors
        'region': 'central'
    })
    
    return df

def create_sample_network_data():
    """Create sample network topology data"""
    # Simple 5-bus network
    buses = pd.DataFrame({
        'bus_id': ['bus_1', 'bus_2', 'bus_3', 'bus_4', 'bus_5'],
        'x': [0, 1, 2, 1, 0.5],
        'y': [0, 0, 0, 1, 0.5],
        'v_nom': [400, 400, 400, 400, 230],  # Voltage levels
        'carrier': ['AC', 'AC', 'AC', 'AC', 'AC']
    })
    
    lines = pd.DataFrame({
        'line_id': ['line_1_2', 'line_2_3', 'line_3_4', 'line_4_5', 'line_5_1'],
        'bus0': ['bus_1', 'bus_2', 'bus_3', 'bus_4', 'bus_5'],
        'bus1': ['bus_2', 'bus_3', 'bus_4', 'bus_5', 'bus_1'],
        'length': [100, 150, 120, 80, 90],  # km
        's_nom': [1000, 800, 600, 400, 500],  # MVA
        'r': [0.01, 0.015, 0.012, 0.008, 0.009],  # Resistance per unit
        'x': [0.1, 0.15, 0.12, 0.08, 0.09]   # Reactance per unit
    })
    
    generators = pd.DataFrame({
        'generator_id': ['coal_1', 'gas_1', 'gas_2', 'nuclear_1', 'hydro_1'],
        'bus': ['bus_1', 'bus_2', 'bus_3', 'bus_1', 'bus_4'],
        'carrier': ['coal', 'gas', 'gas', 'nuclear', 'hydro'],
        'p_nom': [1500, 800, 600, 1200, 400],  # MW
        'marginal_cost': [40, 80, 85, 10, 5],  # $/MWh
        'efficiency': [0.38, 0.55, 0.50, 0.33, 0.90],
        'co2_emissions': [0.82, 0.35, 0.35, 0, 0]  # tons CO2/MWh
    })
    
    return {
        'buses': buses,
        'lines': lines,
        'generators': generators
    }

def main():
    # Create output directories
    os.makedirs('data/raw', exist_ok=True)
    
    # Check if we should use cached data or download fresh data
    use_cached = snakemake.params.data_source.get('use_cached', True)
    
    if use_cached or True:  # For demo, always create sample data
        print("Creating sample data for demonstration...")
        
        # Create sample demand data
        demand_data = create_sample_demand_data()
        demand_data.to_csv(snakemake.output.demand, index=False)
        print(f"Created sample demand data: {snakemake.output.demand}")
        
        # Create sample renewable data
        renewable_data = create_sample_renewable_data()
        renewable_data.to_csv(snakemake.output.renewables, index=False)
        print(f"Created sample renewable data: {snakemake.output.renewables}")
        
        # Create sample network data
        network_data = create_sample_network_data()
        
        # Save network components separately or combine
        combined_network = pd.concat([
            network_data['buses'].assign(component='bus'),
            network_data['lines'].assign(component='line'),
            network_data['generators'].assign(component='generator')
        ], ignore_index=True, sort=False)
        
        combined_network.to_csv(snakemake.output.network, index=False)
        print(f"Created sample network data: {snakemake.output.network}")
        
    else:
        print("Downloading data from external sources...")
        # This would be implemented for real data sources
        # For now, fall back to sample data
        pass
    
    print("Data preparation completed successfully!")

if __name__ == "__main__":
    main()