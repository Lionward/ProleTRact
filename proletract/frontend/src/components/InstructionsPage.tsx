import React, { useMemo, memo } from 'react';
import './InstructionsPage.css';
import './StatisticsDashboard.css';
import './RegionVisualization.css';
import './PopulationComparison.css';

// Color palette for charts - exact same as StatisticsDashboard
const COLOR_PALETTE = [
  '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
  '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6',
  '#60a5fa', '#34d399', '#fbbf24', '#fb7185', '#a78bfa', '#38bdf8',
  '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4', '#84cc16'
];

// Exact replicas of StatisticsDashboard chart components
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
  data: { [key: string]: number };
}

const MotifLengthHistogram: React.FC<MotifLengthHistogramProps> = memo(({ data }) => {
  const entries = Object.entries(data);
  const maxCount = Math.max(...entries.map(([, count]) => count));
  
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

// Exact replica of SequenceTrack from RegionVisualization
interface SequenceTrackProps {
  sequence: string;
  motifIds: string[];
  spans: string;
  label: string;
  metadata: { length: number; motifCoverage: number; motifCount: number; supportingReads?: number };
  motifColors: { [key: number]: string };
  motifNames: string[];
  isReference?: boolean;
}

const SequenceTrack: React.FC<SequenceTrackProps> = ({ 
  sequence, 
  motifIds, 
  spans, 
  label, 
  metadata,
  motifColors,
  motifNames,
  isReference = false
}) => {
  const parseMotifRange = (spans: string): Array<[number, number]> => {
    if (!spans) return [];
    const pattern = /\((\d+)-(\d+)\)/g;
    const ranges: Array<[number, number]> = [];
    let match;
    while ((match = pattern.exec(spans)) !== null) {
      ranges.push([parseInt(match[1]) - 1, parseInt(match[2]) - 1]);
    }
    return ranges;
  };

  const ranges = parseMotifRange(spans);
  const segments: Array<{ start: number; end: number; type: 'motif' | 'interruption'; motifId?: number }> = [];
  
  let pointer = 0;
  ranges.forEach(([start, end], idx) => {
    if (start > pointer) {
      segments.push({ start: pointer, end: start - 1, type: 'interruption' });
    }
    const motifId = idx < motifIds.length ? parseInt(motifIds[idx]) : 0;
    segments.push({ start, end, type: 'motif', motifId });
    pointer = end + 1;
  });
  
  if (pointer < sequence.length) {
    segments.push({ start: pointer, end: sequence.length - 1, type: 'interruption' });
  }

  const charWidth = 8.4;
  const sequenceWidth = sequence.length * charWidth;

  return (
    <div className={`sequence-track ${isReference ? 'reference-track' : ''}`}>
      <div className="sequence-header">
        <div className="sequence-label-container">
          <span className="sequence-label">{label}</span>
        </div>
        <div className="sequence-metadata">
          <div className="metadata-item">
            <span className="metadata-label">Length:</span>
            <span className="metadata-value">{metadata.length} bp</span>
          </div>
          <div className="metadata-item">
            <span className="metadata-label">Motif Coverage:</span>
            <span className="metadata-value">{metadata.motifCoverage}%</span>
          </div>
          <div className="metadata-item">
            <span className="metadata-label">Motif Count:</span>
            <span className="metadata-value">{metadata.motifCount}</span>
          </div>
          {metadata.supportingReads !== undefined && (
            <div className="metadata-item">
              <span className="metadata-label">Reads:</span>
              <span className="metadata-value">{metadata.supportingReads}</span>
            </div>
          )}
        </div>
      </div>
      
      <div className="sequence-scroll-container">
        <div 
          className="sequence-content-wrapper"
          style={{ minWidth: `${Math.max(sequenceWidth, 100)}px` }}
        >
          <div className="sequence-line">
            {segments.map((seg, idx) => {
              const seq = sequence.substring(seg.start, seg.end + 1);
              if (seg.type === 'interruption') {
                return (
                  <span 
                    key={idx} 
                    className="sequence-segment interruption"
                  >
                    {seq.split('').map((base, i) => (
                      <span key={i} className="sequence-base">{base}</span>
                    ))}
                  </span>
                );
              } else {
                const motifId = seg.motifId || 0;
                const color = motifColors[motifId] || '#9ca3af';
                return (
                  <span 
                    key={idx} 
                    className="sequence-segment motif"
                    style={{ backgroundColor: color }}
                  >
                    {seq.split('').map((base, i) => (
                      <span key={i} className="sequence-base">{base}</span>
                    ))}
                  </span>
                );
              }
            })}
          </div>
        </div>
      </div>
    </div>
  );
};

// Exact replica of GenotypeComparisonMatrix from PopulationComparison
interface GenotypeComparisonMatrixProps {
  genotypes: { [sample: string]: string };
}

const GenotypeComparisonMatrix: React.FC<GenotypeComparisonMatrixProps> = memo(({ genotypes }) => {
  const genotypeGroups: { [gt: string]: string[] } = {};
  for (const [sample, gt] of Object.entries(genotypes)) {
    if (!genotypeGroups[gt]) {
      genotypeGroups[gt] = [];
    }
    genotypeGroups[gt].push(sample);
  }
  
  const sortedGroups = Object.entries(genotypeGroups).sort(([a], [b]) => a.localeCompare(b));
  
  const getGenotypeInfo = (gt: string) => {
    if (gt === '0/0') return { icon: '⚪', description: 'Homozygous Reference', color: '#10b981', bgColor: '#d1fae5' };
    if (gt === '0/1') return { icon: '🟡', description: 'Heterozygous', color: '#f59e0b', bgColor: '#fef3c7' };
    if (gt === '1/1') return { icon: '🟠', description: 'Homozygous Alternate', color: '#f97316', bgColor: '#ffedd5' };
    if (gt === '1/2') return { icon: '🔴', description: 'Heterozygous Alternate', color: '#ef4444', bgColor: '#fee2e2' };
    return { icon: '❓', description: 'Unknown', color: '#6b7280', bgColor: '#f3f4f6' };
  };
  
  const totalSamples = Object.keys(genotypes).length;
  const uniqueGenotypes = Object.keys(genotypeGroups).length;
  const homozygous = Object.keys(genotypeGroups).filter(gt => gt === '0/0' || gt === '1/1').length;
  const heterozygous = Object.keys(genotypeGroups).filter(gt => gt === '0/1' || gt === '1/2').length;
  const missing = 0;
  
  return (
    <div className="genotype-matrix">
      <div className="genotype-summary">
        <div className="summary-item">
          <span className="summary-label">Total Samples</span>
          <span className="summary-value">{totalSamples}</span>
        </div>
        <div className="summary-item">
          <span className="summary-label">Unique Genotypes</span>
          <span className="summary-value">{uniqueGenotypes}</span>
        </div>
        <div className="summary-item">
          <span className="summary-label">Homozygous</span>
          <span className="summary-value">{homozygous}</span>
        </div>
        <div className="summary-item">
          <span className="summary-label">Heterozygous</span>
          <span className="summary-value">{heterozygous}</span>
        </div>
        <div className="summary-item">
          <span className="summary-label">Missing</span>
          <span className="summary-value">{missing}</span>
        </div>
      </div>
      
      <div className="genotype-groups">
        {sortedGroups.map(([gt, samples]) => {
          const info = getGenotypeInfo(gt);
          return (
            <div key={gt} className="genotype-group">
              <div className="genotype-group-header">
                <span className="genotype-icon">{info.icon}</span>
                <span className="genotype-value">{gt}</span>
                <span className="genotype-description">- {info.description}</span>
                <span className="genotype-count">({samples.length} samples)</span>
              </div>
              <div className="genotype-samples-tags">
                {samples.map(sample => (
                  <span 
                    key={sample} 
                    className="genotype-tag"
                    data-genotype={gt}
                    style={{ 
                      borderColor: info.color,
                      backgroundColor: info.bgColor,
                      color: info.color
                    }}
                    title={`${sample} - ${gt} (${info.description})`}
                  >
                    {sample}
                  </span>
                ))}
              </div>
            </div>
          );
        })}
      </div>
    </div>
  );
});
GenotypeComparisonMatrix.displayName = 'GenotypeComparisonMatrix';

// Combined Stack Plot and Heatmap for example (exact replica structure)
const ExampleCombinedStackHeatmap: React.FC = () => {
  const samples = ['Sample1', 'Sample2', 'Sample3', 'Sample4'];
  const motifs = ['CAG', 'CAA'];
  const sortedMotifs = ['CAG', 'CAA']; // Sorted by count
  const motifColorMap: { [key: string]: string } = {
    'CAG': COLOR_PALETTE[0],
    'CAA': COLOR_PALETTE[1],
    'Interruption': '#ef4444'
  };
  
  const sampleData = [
    { sample: 'Sample1', segments: [
      { motif: 'CAG', length: 40, start: 0 },
      { motif: 'CAA', length: 3, start: 40 },
      { motif: 'CAG', length: 20, start: 43 }
    ], heatmap: { 'CAG': 13, 'CAA': 1 }},
    { sample: 'Sample2', segments: [
      { motif: 'CAG', length: 50, start: 0 },
      { motif: 'CAG', length: 15, start: 50 }
    ], heatmap: { 'CAG': 21, 'CAA': 0 }},
    { sample: 'Sample3', segments: [
      { motif: 'CAG', length: 30, start: 0 },
      { motif: 'CAA', length: 6, start: 30 },
      { motif: 'CAG', length: 25, start: 36 }
    ], heatmap: { 'CAG': 18, 'CAA': 2 }},
    { sample: 'Sample4', segments: [
      { motif: 'CAG', length: 45, start: 0 }
    ], heatmap: { 'CAG': 15, 'CAA': 0 }}
  ];
  
  const maxLength = 65;
  const maxCount = 21;

  // Heatmap color function (same as CombinedStackHeatmap)
  const getHeatmapColor = (count: number) => {
    const intensity = count / maxCount;
    if (intensity < 0.33) {
      const t = intensity / 0.33;
      const r = Math.round(255 - t * 100);
      const g = Math.round(200 - t * 150);
      const b = Math.round(220 - t * 50);
      return `rgb(${r}, ${g}, ${b})`;
    } else if (intensity < 0.66) {
      const t = (intensity - 0.33) / 0.33;
      const r = Math.round(155 - t * 80);
      const g = Math.round(50 + t * 30);
      const b = Math.round(170 + t * 50);
      return `rgb(${r}, ${g}, ${b})`;
    } else {
      const t = (intensity - 0.66) / 0.34;
      const r = Math.round(75 - t * 30);
      const g = Math.round(80 + t * 20);
      const b = Math.round(220 + t * 35);
      return `rgb(${r}, ${g}, ${b})`;
    }
  };

  return (
    <div className="combined-stack-heatmap">
      {/* Motif Legend at the top */}
      <div className="combined-legend">
        <div className="legend-section">
          <div className="legend-title">Motifs:</div>
          <div className="legend-items-scrollable">
            {motifs.map(motif => (
              <div key={motif} className="legend-item-compact">
                <span className="legend-color-compact" style={{ backgroundColor: motifColorMap[motif] }} />
                <span className="legend-text-compact" title={motif}>
                  {motif}
                </span>
              </div>
            ))}
            <div className="legend-item-compact">
              <span className="legend-color-compact" style={{ backgroundColor: '#ef4444' }} />
              <span className="legend-text-compact">Interruption</span>
            </div>
          </div>
        </div>
      </div>

      <div className="combined-header">
        <div className="combined-sample-label-header">Sample</div>
        <div className="combined-heatmap-header">
          {/* Heatmap Count Legend aligned with heatmap */}
          <div className="heatmap-count-legend">
            <div className="legend-title">Heatmap Count:</div>
            <div className="legend-gradient">
              <span>0</span>
              <div className="gradient-bar" />
              <span>{maxCount}</span>
            </div>
          </div>
        </div>
        <div className="combined-stack-header">Sequence Length</div>
      </div>
      <div className="combined-body">
        {sampleData.map(({ sample, segments, heatmap }) => {
          const totalLength = segments.reduce((sum, s) => sum + s.length, 0);
          return (
            <div key={sample} className="combined-row">
              <div className="combined-sample-label">{sample}</div>
              <div className="combined-heatmap-cells-wrapper">
                <div className="combined-heatmap-cells">
                  {sortedMotifs.map(motif => {
                    const count = (heatmap as { [key: string]: number })[motif] || 0;
                    const color = getHeatmapColor(count);
                    return (
                      <div
                        key={`${sample}_${motif}`}
                        className="heatmap-cell"
                        style={{
                          backgroundColor: count > 0 ? color : '#f3f4f6',
                          opacity: count > 0 ? 1 : 0.6
                        }}
                        title={`${sample} - ${motif}: ${count}`}
                      >
                      </div>
                    );
                  })}
                </div>
              </div>
              <div className="combined-stack-bars-container">
                <div className="combined-stack-bars" style={{ width: `${Math.max((totalLength / maxLength) * 100, 0)}%` }}>
                  {segments.map((segment, idx) => {
                    const widthPercent = totalLength > 0 ? (segment.length / totalLength) * 100 : 0;
                    const color = motifColorMap[segment.motif] || '#9ca3af';
                    return (
                      <div
                        key={idx}
                        className="stack-plot-segment"
                        style={{
                          width: `${widthPercent}%`,
                          backgroundColor: color,
                          borderRight: idx < segments.length - 1 ? '1px solid rgba(0,0,0,0.1)' : 'none'
                        }}
                        title={`${segment.motif}: ${segment.start}-${segment.start + segment.length} (${segment.length}bp)`}
                      />
                    );
                  })}
                </div>
              </div>
              <div className="combined-length-label">{totalLength > 0 ? `${totalLength}bp` : '0bp'}</div>
            </div>
          );
        })}
      </div>
      
      {/* Motif headers at the bottom */}
      <div className="combined-footer">
        <div className="combined-sample-label-footer"></div>
        <div className="combined-heatmap-footer-wrapper">
          <div className="combined-heatmap-footer">
            {sortedMotifs.map(motif => (
              <div key={motif} className="heatmap-motif-footer" title={motif}>
                <span>{motif.length > 6 ? motif.substring(0, 6) + '...' : motif}</span>
              </div>
            ))}
          </div>
        </div>
        <div className="combined-stack-footer"></div>
      </div>
    </div>
  );
};

const InstructionsPage: React.FC = () => {
  // Fake data matching exact structure
  const motifSizeData = { '2-5': 45, '6-10': 120, '11-15': 85, '16-20': 60, '21+': 30 };
  const motifSizeOrder = ['2-5', '6-10', '11-15', '16-20', '21+'];
  const motifLengthData = { '2': 150, '3': 200, '4': 180, '5': 120, '6': 90, '7': 60, '8': 40 };
  const chromosomeData = { 'chr1': 45, 'chr2': 38, 'chr3': 42, 'chr4': 35, 'chr5': 40 };
  const chromosomeOrder = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5'];
  const genotypeData = { '0/0': 180, '0/1': 95, '1/1': 45 };
  const sortedGenotypes: Array<[string, number]> = [['0/0', 180], ['0/1', 95], ['1/1', 45]];

  // Fake sequence data
  const fakeMotifs = ['CAG', 'CAA'];
  const fakeMotifColors: { [key: number]: string } = {
    0: COLOR_PALETTE[0],
    1: COLOR_PALETTE[1]
  };
  const refSequence = 'ATGCAGATGCAGATGCAGATGCAG';
  const refSpans = '(1-3)_(4-6)_(7-9)_(10-12)';
  const refMotifIds = ['0', '0', '0', '0'];
  const altSequence = 'CAGCAGCAGCAGCAGCAGCAGCAA';
  const altSpans = '(1-3)_(4-6)_(7-9)_(10-12)_(13-15)_(16-18)_(19-21)_(22-24)_(25-27)';
  const altMotifIds = ['0', '0', '0', '0', '0', '0', '0', '0', '1'];

  // Fake genotype comparison data
  const fakeGenotypes = {
    'Ref': '0/0',
    'Sample1': '0/1',
    'Sample2': '0/0',
    'Sample3': '1/1',
    'Sample4': '0/1',
    'Sample5': '0/0'
  };

  return (
    <div className="instructions-page">
      <div className="instructions-header">
        <h1>📚 ProleTRact Instructions</h1>
        <p className="instructions-subtitle">
          Learn how to use ProleTRact's visualization and analysis features
        </p>
      </div>

      <div className="instructions-content">
        {/* Statistics Dashboard Section */}
        <section className="instruction-section">
          <h2>📊 Statistics Dashboard</h2>
          <p className="section-description">
            The Statistics Dashboard provides comprehensive overview of your VCF file data.
          </p>
          
          <div className="plot-example">
            <h3>Available Charts:</h3>
            <div className="example-plots-grid">
              <div className="chart-container">
                <h3>Regions by Max Motif Size</h3>
                <MotifSizeBarChart data={motifSizeData} categoryOrder={motifSizeOrder} />
              </div>
              <div className="chart-container">
                <h3>Distribution of Motif Sizes (All Motifs)</h3>
                <MotifLengthHistogram data={motifLengthData} />
              </div>
              <div className="chart-container full-width">
                <h3>Regions by Chromosome</h3>
                <ChromosomeBarChart data={chromosomeData} sortedChromosomes={chromosomeOrder} />
              </div>
              <div className="chart-container full-width">
                <h3>Genotype Distribution</h3>
                <GenotypeBarChart data={genotypeData} sortedGenotypes={sortedGenotypes} />
              </div>
            </div>
            <div className="usage-tip">
              <strong>💡 Tip:</strong> Access the Statistics Dashboard by clicking the "Statistics" tab 
              after loading a VCF file.
            </div>
          </div>
        </section>

        {/* Region Visualization Section */}
        <section className="instruction-section">
          <h2>🔬 Region Visualization</h2>
          <p className="section-description">
            Detailed sequence-level visualization of tandem repeat regions with color-coded motifs.
          </p>
          
          <div className="plot-example">
            <h3>Features:</h3>
            <div className="example-sequence-container">
              <div className="motif-legend">
                <div className="legend-header">Motif Legend</div>
                <div className="legend-items">
                  {fakeMotifs.map((motif, idx) => (
                    <div key={idx} className="legend-item">
                      <span
                        className="legend-color"
                        style={{ backgroundColor: fakeMotifColors[idx] }}
                      />
                      <span className="legend-text">{motif} ({motif.length}bp)</span>
                    </div>
                  ))}
                </div>
              </div>
              <SequenceTrack
                sequence={refSequence}
                motifIds={refMotifIds}
                spans={refSpans}
                label="Reference"
                metadata={{ length: 24, motifCoverage: 50, motifCount: 4 }}
                motifColors={fakeMotifColors}
                motifNames={fakeMotifs}
                isReference={true}
              />
              <SequenceTrack
                sequence={altSequence}
                motifIds={altMotifIds}
                spans={altSpans}
                label="Alternate Allele"
                metadata={{ length: 27, motifCoverage: 100, motifCount: 9, supportingReads: 42 }}
                motifColors={fakeMotifColors}
                motifNames={fakeMotifs}
              />
            </div>
            <ul className="feature-list">
              <li>
                <strong>Sequence Display:</strong> Visual representation of reference and alternate alleles 
                with highlighted motif regions. Each motif is color-coded for easy identification.
              </li>
              <li>
                <strong>Motif Legend:</strong> Color legend showing which color corresponds to which motif. 
                Hover over sequence segments to see detailed information.
              </li>
              <li>
                <strong>Metadata Display:</strong> Shows region length, motif coverage, motif count, and 
                supporting reads information.
              </li>
              <li>
                <strong>Horizontal Scrolling:</strong> Long sequences can be scrolled horizontally to view 
                the entire region.
              </li>
              <li>
                <strong>Scale Bars:</strong> Visual scale indicators showing sequence length and position.
              </li>
            </ul>
            <div className="usage-tip">
              <strong>💡 Tip:</strong> Select a region from the sidebar navigation panel to view its detailed 
              visualization. Use the search bar to quickly find specific regions.
            </div>
          </div>
        </section>

        {/* Pathogenic Catalog Section */}
        <section className="instruction-section">
          <h2>⚠️ Pathogenic Catalog Integration</h2>
          <p className="section-description">
            Identify and filter regions with pathogenic tandem repeats using integrated pathogenic catalog data.
          </p>
          
          <div className="plot-example">
            <h3>Features:</h3>
            <ul className="feature-list">
              <li>
                <strong>Pathogenic Threshold Display:</strong> When a region has a defined pathogenic threshold, 
                a red dashed line appears on bar plots and stack plots indicating the threshold value. The Y-axis 
                automatically extends to include the threshold if it exceeds the maximum data value.
              </li>
              <li>
                <strong>Pathogenic Filter:</strong> Enable the "Show only regions with pathogenic alleles" checkbox 
                to filter the region list to only show regions where at least one allele exceeds the pathogenic threshold. 
                When active, navigation, search, and all region listings only show pathogenic regions.
              </li>
              <li>
                <strong>Gene Search:</strong> Search for regions by gene name in the search bar. The system will 
                automatically find regions associated with the gene from the pathogenic catalog.
              </li>
              <li>
                <strong>Pathogenicity Panel:</strong> When viewing a region with pathogenic catalog information, 
                a panel displays the gene name, associated disease, inheritance pattern, and pathogenic threshold. 
                The panel is color-coded: red if any allele exceeds the threshold, green if below threshold, 
                or neutral if no threshold is defined.
              </li>
              <li>
                <strong>Allele Color Coding:</strong> In individual mode, alleles are color-coded based on 
                pathogenicity:
                <ul className="nested-list">
                  <li><strong>Red/Coral border:</strong> Pathogenic allele (exceeds threshold)</li>
                  <li><strong>Green border:</strong> Non-pathogenic allele (below threshold)</li>
                  <li><strong>No border:</strong> Reference allele or no threshold defined</li>
                </ul>
              </li>
              <li>
                <strong>Cohort Mode Highlighting:</strong> In cohort analysis, samples that exceed the pathogenic 
                threshold are highlighted with a red border and shadow around the entire sample row in the stack plot.
              </li>
            </ul>
            <div className="usage-tip">
              <strong>💡 Tip:</strong> The pathogenic catalog is automatically loaded from a BED file. Regions with 
              pathogenic information will show the threshold line on plots and allow filtering. Use gene search to 
              quickly find regions of clinical interest.
            </div>
          </div>
        </section>

        {/* Population Comparison Section */}
        <section className="instruction-section">
          <h2>🌐 Population Comparison</h2>
          <p className="section-description">
            Compare your sample against population data to identify rare or novel expansions.
          </p>
          
          <div className="plot-example">
            <h3>Available Visualizations:</h3>
            <GenotypeComparisonMatrix genotypes={fakeGenotypes} />
            <ul className="feature-list">
              <li>
                <strong>Genotype Comparison Matrix:</strong> Grid view showing genotype comparisons between 
                your sample and population samples. Color-coded to highlight differences.
              </li>
              <li>
                <strong>Sequence Comparison:</strong> Side-by-side comparison of sequences from your sample 
                and population samples, showing motif patterns and differences.
              </li>
              <li>
                <strong>Genotype Summary Statistics:</strong> Summary cards showing counts of homozygous, 
                heterozygous, and missing genotypes in the population.
              </li>
              <li>
                <strong>Pathogenic Threshold Display:</strong> Bar plots show a red dashed line indicating 
                the pathogenic threshold when available for the region.
              </li>
            </ul>
            <div className="usage-tip">
              <strong>💡 Tip:</strong> Load a public VCF folder containing population data from the sidebar 
              to enable population comparison features.
            </div>
          </div>
        </section>

        {/* Cohort Analysis Section */}
        <section className="instruction-section">
          <h2>👥 Cohort Analysis</h2>
          <p className="section-description">
            Analyze multiple samples simultaneously to identify patterns and variations across a cohort.
          </p>
          
          <div className="plot-example">
            <h3>Available Plots:</h3>
            <div className="example-plots-grid">
              <div className="chart-container full-width">
                <h3>Genotype Frequency Distribution</h3>
                <GenotypeBarChart data={genotypeData} sortedGenotypes={sortedGenotypes} />
              </div>
              <div className="chart-container full-width">
                <h3>Motif Stack Plot & Heatmap</h3>
                <ExampleCombinedStackHeatmap />
              </div>
            </div>
            <ul className="feature-list">
              <li>
                <strong>Genotype Frequency Distribution:</strong> Bar chart showing the frequency of each 
                genotype across all samples in the cohort. Helps identify common vs. rare genotypes.
              </li>
              <li>
                <strong>Motif Stack Plot & Heatmap:</strong> 
                <ul className="nested-list">
                  <li><strong>Stack Plot:</strong> Visual representation of motif runs per sample/allele, 
                  stacked to show cumulative patterns.</li>
                  <li><strong>Heatmap:</strong> Color-coded matrix showing motif occurrence counts across 
                  samples. Darker colors indicate higher counts.</li>
                </ul>
              </li>
              <li>
                <strong>Bar Plot for Motif Count:</strong> Bar chart showing total motif segments per 
                sample/allele. Includes pathogenic threshold lines (red dashed line) when available for 
                clinical interpretation. The Y-axis automatically extends to include the threshold value.
              </li>
              <li>
                <strong>Cluster Plot:</strong> K-means clustering visualization grouping samples by 
                similarity in copy number, length, or both. Helps identify sample groups with similar 
                repeat characteristics.
              </li>
            </ul>
            <div className="usage-tip">
              <strong>💡 Tip:</strong> Switch to "Cohort Read" or "Cohort Assembly" mode from the sidebar 
              to access cohort analysis features. Load a folder containing multiple VCF files.
            </div>
          </div>
        </section>

        {/* Getting Started Section */}
        <section className="instruction-section">
          <h2>🚀 Getting Started</h2>
          <div className="steps-list">
            <div className="step-item">
              <div className="step-number">1</div>
              <div className="step-content">
                <h3>Load Your VCF File</h3>
                <p>Enter the path to your TandemTwister-generated VCF file in the sidebar and click "Load VCF". 
                The file should be bgzipped and tabix-indexed.</p>
              </div>
            </div>
            <div className="step-item">
              <div className="step-number">2</div>
              <div className="step-content">
                <h3>Filter by Genotype (Optional)</h3>
                <p>Use the genotype filter in the sidebar to focus on specific genotypes of interest.</p>
              </div>
            </div>
            <div className="step-item">
              <div className="step-number">3</div>
              <div className="step-content">
                <h3>Select a Region</h3>
                <p>Choose a region from the navigation panel, use the search functionality to find specific regions, 
                or search by gene name to find regions associated with a particular gene.</p>
              </div>
            </div>
            <div className="step-item">
              <div className="step-number">4</div>
              <div className="step-content">
                <h3>Use Pathogenic Filter (Optional)</h3>
                <p>Enable the "Show only regions with pathogenic alleles" checkbox to filter to only regions where 
                at least one allele exceeds the pathogenic threshold. This affects navigation, search, and all region listings.</p>
              </div>
            </div>
            <div className="step-item">
              <div className="step-number">5</div>
              <div className="step-content">
                <h3>Explore Visualizations</h3>
                <p>Switch between the Statistics and Region Analysis tabs to explore different views of your data. 
                Look for pathogenic threshold lines on plots and check the Pathogenicity Panel for clinical information.</p>
              </div>
            </div>
            <div className="step-item">
              <div className="step-number">6</div>
              <div className="step-content">
                <h3>Compare with Population Data</h3>
                <p>Load a public VCF folder to compare your sample against population data and identify rare variants.</p>
              </div>
            </div>
          </div>
        </section>

        {/* Keyboard Shortcuts Section */}
        <section className="instruction-section">
          <h2>⌨️ Keyboard Shortcuts</h2>
          <div className="shortcuts-grid">
            <div className="shortcut-item">
              <kbd>/</kbd>
              <span>Focus search bar</span>
            </div>
            <div className="shortcut-item">
              <kbd>?</kbd>
              <span>Show keyboard shortcuts help</span>
            </div>
          </div>
        </section>
      </div>
    </div>
  );
};

export default InstructionsPage;
