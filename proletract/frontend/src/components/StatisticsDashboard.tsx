//implement the following: an option to add a folder where the population analysis is in the indiviudal mode, for the tab allele vs population where proleTRact shows normally genotype information, heatmap and stackplot. cluster plot and bar plot
import React, { useEffect, useState, useMemo, memo, useRef } from 'react';
import axios from 'axios';
import './StatisticsDashboard.css';
import ExportMenu from './ExportMenu';
import { exportToJSON, generateFilename } from '../utils/exportUtils';

// Colorful palette for charts
const COLOR_PALETTE = [
  '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
  '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6',
  '#60a5fa', '#34d399', '#fbbf24', '#fb7185', '#a78bfa', '#38bdf8',
  '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4', '#84cc16'
];

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

interface StatisticsDashboardProps {
  vcfPath: string;
}

interface VCFStatistics {
  total_regions: number;
  num_chromosomes: number;
  avg_motif_size: number;
  max_motif_size: number;
  motif_size_counts: { [key: string]: number };
  motif_length_histogram: { [key: string]: number };  
  regions_by_chromosome: { [key: string]: number };
  genotype_counts: { [key: string]: number };
}

const StatisticsDashboard: React.FC<StatisticsDashboardProps> = ({ vcfPath }) => {
  const [stats, setStats] = useState<VCFStatistics | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [cachedVcfPath, setCachedVcfPath] = useState<string>('');

  useEffect(() => {
    if (!vcfPath) {
      setStats(null);
      setCachedVcfPath('');
      return;
    }

    // Only fetch if VCF path changed (component keeps data when just hidden)
    if (vcfPath === cachedVcfPath && stats !== null) {
      return; // Already have data for this VCF file, no need to refetch
    }

    const fetchStatistics = async () => {
      setLoading(true);
      setError(null);
      try {
        const response = await axios.get(`${API_BASE}/api/vcf/statistics`, {
          params: { vcf_path: vcfPath }
        });
        setStats(response.data);
        setCachedVcfPath(vcfPath);
      } catch (err: any) {
        setError(err.response?.data?.detail || 'Failed to load statistics');
        setStats(null);
        setCachedVcfPath('');
      } finally {
        setLoading(false);
      }
    };

    fetchStatistics();
  }, [vcfPath]); // Only depend on vcfPath - component state persists when hidden

  if (loading) {
    return <div className="loading-message">Loading statistics...</div>;
  }

  if (error) {
    return <div className="error-message">Error: {error}</div>;
  }

  if (!stats) {
    return <div className="no-data-message">No statistics available</div>;
  }

  return <StatisticsDisplay stats={stats} />;
};

interface StatisticsDisplayProps {
  stats: VCFStatistics;
}

const StatisticsDisplay: React.FC<StatisticsDisplayProps> = memo(({ stats }) => {
  const dashboardRef = useRef<HTMLDivElement>(null);
  
  // Sort chromosomes naturally
  const sortedChromosomes = Object.keys(stats.regions_by_chromosome).sort((a, b) => {
    const aNum = a.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25');
    const bNum = b.replace('chr', '').replace('X', '23').replace('Y', '24').replace('M', '25');
    const aInt = parseInt(aNum) || 999;
    const bInt = parseInt(bNum) || 999;
    return aInt - bInt;
  });

  // Sort genotype counts
  const sortedGenotypes = Object.entries(stats.genotype_counts).sort((a, b) => {
    const aParts = a[0].split('/').length;
    const bParts = b[0].split('/').length;
    if (aParts !== bParts) return aParts - bParts;
    return a[0].localeCompare(b[0]);
  });

  // Prepare motif size categories
  const categoryOrder = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", ">10", "Unknown"];
  const sortedMotifSizes = categoryOrder.filter(cat => stats.motif_size_counts[cat] !== undefined);

  const handleExportJSON = () => {
    exportToJSON(stats, generateFilename('vcf_statistics', 'json'));
  };

  return (
    <div className="statistics-dashboard" ref={dashboardRef}>
      <div className="statistics-header">
        <h2>📊 VCF Statistics</h2>
        <ExportMenu
          elementRef={dashboardRef}
          filename="vcf_statistics"
          onExportJSON={handleExportJSON}
          showDataExport={true}
        />
      </div>

      {/* Summary Metrics */}
      <div className="summary-metrics">
        <div className="metric-card">
          <div className="metric-label">Total Regions</div>
          <div className="metric-value">{stats.total_regions.toLocaleString()}</div>
        </div>
        <div className="metric-card">
          <div className="metric-label">Chromosomes</div>
          <div className="metric-value">{stats.num_chromosomes}</div>
        </div>
        <div className="metric-card">
          <div className="metric-label">Avg Motif Size</div>
          <div className="metric-value">{stats.avg_motif_size.toFixed(1)} bp</div>
        </div>
        <div className="metric-card">
          <div className="metric-label">Max Motif Size</div>
          <div className="metric-value">{stats.max_motif_size} bp</div>
        </div>
      </div>

      {/* Genotype Distribution */}
      {sortedGenotypes.length > 0 && (
        <div className="statistics-section">
          <h3>🧬 Genotype Distribution</h3>
          <div className="genotype-metrics">
            {sortedGenotypes.slice(0, 6).map(([gt, count]) => {
              const percentage = ((count / stats.total_regions) * 100).toFixed(1);
              return (
                <div key={gt} className="genotype-card">
                  <div className="genotype-label">Genotype {gt}</div>
                  <div className="genotype-value">{count.toLocaleString()}</div>
                  <div className="genotype-percentage">{percentage}%</div>
                </div>
              );
            })}
          </div>
          {sortedGenotypes.length > 6 && (
            <div className="other-genotypes">
              <h4>Other Genotypes:</h4>
              <div className="genotype-table">
                <table>
                  <thead>
                    <tr>
                      <th>Genotype</th>
                      <th>Count</th>
                      <th>Percentage</th>
                    </tr>
                  </thead>
                  <tbody>
                    {sortedGenotypes.slice(6).map(([gt, count]) => {
                      const percentage = ((count / stats.total_regions) * 100).toFixed(2);
                      return (
                        <tr key={gt}>
                          <td>{gt}</td>
                          <td>{count.toLocaleString()}</td>
                          <td>{percentage}%</td>
                        </tr>
                      );
                    })}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      )}

      {/* Charts Section */}
      <div className="charts-section">
        <div className="chart-row">
          {/* Motif Size Distribution */}
          <div className="chart-container">
            <h3>Regions by Max Motif Size</h3>
            <MotifSizeBarChart data={stats.motif_size_counts} categoryOrder={sortedMotifSizes} />
          </div>

          {/* Motif Length Distribution */}
          <div className="chart-container">
            <h3>Distribution of Motif Sizes (All Motifs)</h3>
            <p style={{ fontSize: '0.85rem', color: '#6b7280', marginTop: '-0.5rem', marginBottom: '1rem' }}>
              Shows the count of individual motifs by size. Multiple motifs per region are counted separately.
            </p>
            <MotifLengthHistogram data={stats.motif_length_histogram} />
          </div>
        </div>

        {/* Chromosome Distribution */}
        <div className="chart-container full-width">
          <h3>Regions by Chromosome</h3>
          <ChromosomeBarChart 
            data={stats.regions_by_chromosome} 
            sortedChromosomes={sortedChromosomes}
          />
        </div>

        {/* Genotype Distribution Chart */}
        {sortedGenotypes.length > 0 && (
          <div className="chart-container full-width">
            <h3>Genotype Distribution</h3>
            <GenotypeBarChart data={stats.genotype_counts} sortedGenotypes={sortedGenotypes} />
          </div>
        )}
      </div>
    </div>
  );
});
StatisticsDisplay.displayName = 'StatisticsDisplay';

// Chart Components
interface MotifSizeBarChartProps {
  data: { [key: string]: number };
  categoryOrder: string[];
}

const MotifSizeBarChart: React.FC<MotifSizeBarChartProps> = memo(({ data, categoryOrder }) => {
  const maxValue = Math.max(...Object.values(data));
  
  return (
    <div className="bar-chart vertical">
      <div className="bar-chart-container">
        {categoryOrder.map((category, index) => {
          const value = data[category] || 0;
          const percentage = maxValue > 0 ? (value / maxValue) * 100 : 0;
          const color = COLOR_PALETTE[index % COLOR_PALETTE.length];
          return (
            <div key={category} className="bar-item-vertical">
              <div className="bar-wrapper-vertical">
                <div 
                  className="bar-fill-vertical"
                  style={{ 
                    height: `${percentage}%`,
                    backgroundColor: color
                  }}
                  title={`${category}: ${value.toLocaleString()}`}
                >
                  <span className="bar-value-vertical">{value.toLocaleString()}</span>
                </div>
              </div>
              <div className="bar-label-vertical">{category}</div>
            </div>
          );
        })}
      </div>
    </div>
  );
});
MotifSizeBarChart.displayName = 'MotifSizeBarChart';

interface MotifLengthHistogramProps {
  data: { [key: string]: number };  // Changed to accept binned histogram data
}

const MotifLengthHistogram: React.FC<MotifLengthHistogramProps> = memo(({ data }) => {
  if (!data || Object.keys(data).length === 0) {
    return <div className="no-data">No motif length data available</div>;
  }

  // Data is already binned from backend, just need to sort and display
  const entries = Object.entries(data);
  const maxCount = Math.max(...entries.map(([, count]) => count));
  
  if (maxCount === 0) {
    return <div className="no-data">No data to display</div>;
  }

  // Sort bins by their numeric range (parse the label like "10-20" or "10")
  const sortedEntries = entries.sort(([labelA], [labelB]) => {
    const numA = parseInt(labelA.split('-')[0]) || 0;
    const numB = parseInt(labelB.split('-')[0]) || 0;
    return numA - numB;
  });

  return (
    <div className="histogram-container">
      <div className="histogram-chart">
        {sortedEntries.map(([label, count], index) => {
          const percentage = maxCount > 0 ? (count / maxCount) * 100 : 0;
          // Use a consistent color scheme for histogram (alternating colors for better distinction)
          const colorIndex = index % 12;
          const color = COLOR_PALETTE[colorIndex];
          return (
            <div key={label} className="histogram-bar-item">
              <div className="histogram-bar-wrapper">
                <div 
                  className="histogram-bar"
                  style={{ 
                    height: `${percentage}%`,
                    backgroundColor: color
                  }}
                  title={`${label} bp: ${count.toLocaleString()} motifs`}
                >
                  {percentage > 8 && (
                    <span className="histogram-bar-value">{count.toLocaleString()}</span>
                  )}
                </div>
              </div>
              <div className="histogram-bar-label">{label}</div>
            </div>
          );
        })}
      </div>
    </div>
  );
});
MotifLengthHistogram.displayName = 'MotifLengthHistogram';

interface ChromosomeBarChartProps {
  data: { [key: string]: number };
  sortedChromosomes: string[];
}

const ChromosomeBarChart: React.FC<ChromosomeBarChartProps> = memo(({ data, sortedChromosomes }) => {
  const maxValue = Math.max(...Object.values(data));
  
  return (
    <div className="bar-chart vertical">
      <div className="bar-chart-container">
        {sortedChromosomes.map((chrom, index) => {
          const value = data[chrom] || 0;
          const percentage = maxValue > 0 ? (value / maxValue) * 100 : 0;
          const color = COLOR_PALETTE[index % COLOR_PALETTE.length];
          return (
            <div key={chrom} className="bar-item-vertical">
              <div className="bar-wrapper-vertical">
                <div 
                  className="bar-fill-vertical"
                  style={{ 
                    height: `${percentage}%`,
                    backgroundColor: color
                  }}
                  title={`${chrom}: ${value.toLocaleString()}`}
                >
                  <span className="bar-value-vertical">{value.toLocaleString()}</span>
                </div>
              </div>
              <div className="bar-label-vertical">{chrom}</div>
            </div>
          );
        })}
      </div>
    </div>
  );
});
ChromosomeBarChart.displayName = 'ChromosomeBarChart';

interface GenotypeBarChartProps {
  data: { [key: string]: number };
  sortedGenotypes: Array<[string, number]>;
}

const GenotypeBarChart: React.FC<GenotypeBarChartProps> = memo(({ data, sortedGenotypes }) => {
  const maxValue = Math.max(...Object.values(data));
  
  return (
    <div className="bar-chart vertical">
      <div className="bar-chart-container">
        {sortedGenotypes.map(([gt, count], index) => {
          const percentage = maxValue > 0 ? (count / maxValue) * 100 : 0;
          const color = COLOR_PALETTE[index % COLOR_PALETTE.length];
          return (
            <div key={gt} className="bar-item-vertical">
              <div className="bar-wrapper-vertical">
                <div 
                  className="bar-fill-vertical"
                  style={{ 
                    height: `${percentage}%`,
                    backgroundColor: color
                  }}
                  title={`${gt}: ${count.toLocaleString()}`}
                >
                  <span className="bar-value-vertical">{count.toLocaleString()}</span>
                </div>
              </div>
              <div className="bar-label-vertical">{gt}</div>
            </div>
          );
        })}
      </div>
    </div>
  );
});
GenotypeBarChart.displayName = 'GenotypeBarChart';

export default StatisticsDashboard;

