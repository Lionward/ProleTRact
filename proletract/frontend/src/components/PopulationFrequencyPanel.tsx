import React, { useState, useEffect, useMemo } from 'react';
import axios from 'axios';
import './PopulationFrequencyPanel.css';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

interface AlleleFrequency {
  allele: string;
  count: number;
  frequency: number;
  percentage: number;
}

interface PopulationFrequencyData {
  region: string;
  total_samples: number;
  allele_frequencies: AlleleFrequency[];
  rare_alleles: string[];
  novel_alleles: string[];
}

interface PopulationFrequencyPanelProps {
  chr: string;
  pos: number;
  stop: number;
  region: string;
  currentAlleles?: string[];
}

const PopulationFrequencyPanel: React.FC<PopulationFrequencyPanelProps> = ({
  chr,
  pos,
  stop,
  region,
  currentAlleles = [],
}) => {
  const [frequencyData, setFrequencyData] = useState<PopulationFrequencyData | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const fetchFrequencyData = async () => {
      setLoading(true);
      setError(null);
      try {
        // This endpoint would need to be implemented in the backend
        const response = await axios.get(`${API_BASE}/api/population/frequency`, {
          params: {
            chr: chr,
            start: pos,
            end: stop,
          },
        });
        setFrequencyData(response.data);
      } catch (err: any) {
        // Silently fail if endpoint doesn't exist yet
        if (err.response?.status !== 404) {
          console.error('Error fetching population frequency:', err);
          setError(err.response?.data?.detail || 'Failed to load population frequency data');
        }
        setFrequencyData(null);
      } finally {
        setLoading(false);
      }
    };

    if (chr && pos && stop) {
      fetchFrequencyData();
    }
  }, [chr, pos, stop]);

  const maxFrequency = useMemo(() => {
    if (!frequencyData) return 0;
    return Math.max(...frequencyData.allele_frequencies.map(a => a.frequency), 0);
  }, [frequencyData]);

  if (loading) {
    return (
      <div className="population-frequency-panel loading">
        <div className="population-frequency-header">
          <h3>📊 Population Frequency</h3>
        </div>
        <div className="population-frequency-content">
          <p>Loading population frequency data...</p>
        </div>
      </div>
    );
  }

  if (error || !frequencyData) {
    return null; // Don't show panel if there's an error or no data
  }

  const isRare = (allele: string) => frequencyData.rare_alleles.includes(allele);
  const isNovel = (allele: string) => frequencyData.novel_alleles.includes(allele);
  const isCurrent = (allele: string) => currentAlleles.includes(allele);

  return (
    <div className="population-frequency-panel">
      <div className="population-frequency-header">
        <h3>📊 Population Frequency</h3>
        <div className="population-frequency-summary">
          <span className="frequency-summary-item">
            Total Samples: <strong>{frequencyData.total_samples}</strong>
          </span>
          <span className="frequency-summary-item">
            Alleles: <strong>{frequencyData.allele_frequencies.length}</strong>
          </span>
        </div>
      </div>

      <div className="population-frequency-content">
        <div className="frequency-chart-container">
          <h4>Allele Frequency Distribution</h4>
          <div className="frequency-chart">
            {frequencyData.allele_frequencies
              .sort((a, b) => b.frequency - a.frequency)
              .map((alleleFreq, index) => {
                const widthPercent = maxFrequency > 0 ? (alleleFreq.frequency / maxFrequency) * 100 : 0;
                const isRareAllele = isRare(alleleFreq.allele);
                const isNovelAllele = isNovel(alleleFreq.allele);
                const isCurrentAllele = isCurrent(alleleFreq.allele);

                return (
                  <div
                    key={alleleFreq.allele}
                    className={`frequency-bar-item ${isCurrentAllele ? 'current-allele' : ''} ${isRareAllele ? 'rare' : ''} ${isNovelAllele ? 'novel' : ''}`}
                  >
                    <div className="frequency-bar-label">
                      <span className="allele-label">
                        {alleleFreq.allele}
                        {isCurrentAllele && <span className="current-badge">Current</span>}
                        {isRareAllele && <span className="rare-badge">Rare</span>}
                        {isNovelAllele && <span className="novel-badge">Novel</span>}
                      </span>
                      <span className="frequency-value">
                        {alleleFreq.count} ({alleleFreq.percentage.toFixed(2)}%)
                      </span>
                    </div>
                    <div className="frequency-bar-wrapper">
                      <div
                        className="frequency-bar"
                        style={{
                          width: `${widthPercent}%`,
                          backgroundColor: isCurrentAllele
                            ? '#667eea'
                            : isRareAllele
                            ? '#f59e0b'
                            : isNovelAllele
                            ? '#ef4444'
                            : '#4fd1c7',
                        }}
                        title={`${alleleFreq.allele}: ${alleleFreq.count} samples (${alleleFreq.percentage.toFixed(2)}%)`}
                      />
                    </div>
                  </div>
                );
              })}
          </div>
        </div>

        {frequencyData.rare_alleles.length > 0 && (
          <div className="frequency-alerts">
            <div className="frequency-alert rare-alert">
              <strong>⚠️ Rare Alleles Detected:</strong>
              <div className="alert-alleles">
                {frequencyData.rare_alleles.map(allele => (
                  <span key={allele} className="alert-allele">
                    {allele}
                  </span>
                ))}
              </div>
            </div>
          </div>
        )}

        {frequencyData.novel_alleles.length > 0 && (
          <div className="frequency-alerts">
            <div className="frequency-alert novel-alert">
              <strong>🔬 Novel Alleles Detected:</strong>
              <div className="alert-alleles">
                {frequencyData.novel_alleles.map(allele => (
                  <span key={allele} className="alert-allele">
                    {allele}
                  </span>
                ))}
              </div>
            </div>
          </div>
        )}

        <div className="frequency-legend">
          <div className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#667eea' }} />
            <span>Current Sample</span>
          </div>
          <div className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#f59e0b' }} />
            <span>Rare (&lt;1%)</span>
          </div>
          <div className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#ef4444' }} />
            <span>Novel (Not in database)</span>
          </div>
          <div className="legend-item">
            <span className="legend-color" style={{ backgroundColor: '#4fd1c7' }} />
            <span>Common</span>
          </div>
        </div>
      </div>
    </div>
  );
};

export default PopulationFrequencyPanel;




