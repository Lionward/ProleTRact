import React, { useState, useEffect } from 'react';
import axios from 'axios';
import CohortAnalysis from './CohortAnalysis';
import './RegionVisualization.css';
import './FloatingNavigation.css';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

type AppMode = 'cohort-read' | 'cohort-assembly';

interface CohortModeViewProps {
  mode: AppMode;
  cohortFolder: string;
  publicVcfFolder: string;
  loading: boolean;
  selectedRegion?: string;
  onRegionSelect?: (region: string) => void;
}

const CohortModeView: React.FC<CohortModeViewProps> = ({
  mode,
  cohortFolder,
  publicVcfFolder,
  loading,
  selectedRegion = '',
  onRegionSelect
}) => {
  const [availableRegions, setAvailableRegions] = useState<string[]>([]);
  const [loadingRegions, setLoadingRegions] = useState(false);

  // Determine the actual cohort folder path
  const cohortFolderPath = publicVcfFolder || cohortFolder;

  // Fetch available regions for navigation (but not for search - that's in sidebar)
  useEffect(() => {
    if (!cohortFolderPath) {
      setAvailableRegions([]);
      return;
    }

    const fetchRegions = async () => {
      setLoadingRegions(true);
      try {
        const response = await axios.get(`${API_BASE}/api/population/regions`, {
          params: { folder_path: cohortFolderPath }
        });
        const regions = response.data.regions || [];
        // Sort regions properly by chromosome and position
        const sortedRegions = [...regions].sort((a, b) => {
          const parseRegion = (r: string) => {
            const match = r.match(/^([^:]+):(\d+)-(\d+)$/);
            if (!match) return { chr: '', start: 0, end: 0 };
            return {
              chr: match[1],
              start: parseInt(match[2], 10),
              end: parseInt(match[3], 10)
            };
          };
          const aParsed = parseRegion(a);
          const bParsed = parseRegion(b);
          
          // Compare chromosome
          if (aParsed.chr !== bParsed.chr) {
            // Handle numeric vs non-numeric chromosomes
            const aNum = parseInt(aParsed.chr.replace(/^chr/i, ''), 10);
            const bNum = parseInt(bParsed.chr.replace(/^chr/i, ''), 10);
            if (!isNaN(aNum) && !isNaN(bNum)) {
              return aNum - bNum;
            }
            return aParsed.chr.localeCompare(bParsed.chr);
          }
          
          // Same chromosome, compare start position
          return aParsed.start - bParsed.start;
        });
        setAvailableRegions(sortedRegions);
      } catch (err: any) {
        console.error('Failed to fetch regions:', err);
        setAvailableRegions([]);
      } finally {
        setLoadingRegions(false);
      }
    };

    fetchRegions();
  }, [cohortFolderPath]);

  // When selectedRegion changes, ensure it's in availableRegions list
  useEffect(() => {
    if (!selectedRegion || !cohortFolderPath) return;
    
    // If the region is not in the list, add it and sort
    // This ensures the navigator can show the correct position
    setAvailableRegions(prev => {
      // Check if the selected region is already in the list
      if (prev.includes(selectedRegion)) {
        return prev; // Already in the list
      }
      
      const updated = [...prev, selectedRegion];
      // Sort regions to maintain order (chromosome, then position)
      updated.sort((a, b) => {
        const parseRegion = (r: string) => {
          const match = r.match(/^([^:]+):(\d+)-(\d+)$/);
          if (!match) return { chr: '', start: 0, end: 0 };
          return {
            chr: match[1],
            start: parseInt(match[2], 10),
            end: parseInt(match[3], 10)
          };
        };
        const aParsed = parseRegion(a);
        const bParsed = parseRegion(b);
        
        // Compare chromosome
        if (aParsed.chr !== bParsed.chr) {
          // Handle numeric vs non-numeric chromosomes
          const aNum = parseInt(aParsed.chr.replace(/^chr/i, ''), 10);
          const bNum = parseInt(bParsed.chr.replace(/^chr/i, ''), 10);
          if (!isNaN(aNum) && !isNaN(bNum)) {
            return aNum - bNum;
          }
          return aParsed.chr.localeCompare(bParsed.chr);
        }
        
        // Same chromosome, compare start position
        return aParsed.start - bParsed.start;
      });
      return updated;
    });
  }, [selectedRegion, cohortFolderPath]);

  // Get current region index and navigation functions
  const currentRegionIndex = availableRegions.findIndex(region => region === selectedRegion);
  const hasPreviousRegion = currentRegionIndex > 0;
  const hasNextRegion = currentRegionIndex >= 0 && currentRegionIndex < availableRegions.length - 1;

  const handlePreviousRegion = () => {
    if (hasPreviousRegion && currentRegionIndex > 0) {
      const prevRegion = availableRegions[currentRegionIndex - 1];
      if (onRegionSelect) {
        onRegionSelect(prevRegion);
      }
    }
  };

  const handleNextRegion = () => {
    if (hasNextRegion && currentRegionIndex < availableRegions.length - 1) {
      const nextRegion = availableRegions[currentRegionIndex + 1];
      if (onRegionSelect) {
        onRegionSelect(nextRegion);
      }
    }
  };

  if (!cohortFolderPath) {
    return (
      <div className="welcome-screen">
        <div className="welcome-content">
          <div className="welcome-icon">👥</div>
          <h2>Cohort Mode</h2>
          <p>Analyze tandem repeat regions across multiple samples from a cohort folder</p>
          <div className="welcome-features">
            <div className="feature-card">
              <span className="feature-icon">📁</span>
              <h3>Load Cohort Folder</h3>
              <p>Load a folder containing multiple VCF files (h1.vcf.gz, h2.vcf.gz)</p>
            </div>
            <div className="feature-card">
              <span className="feature-icon">🔍</span>
              <h3>Search Regions</h3>
              <p>Enter a region to analyze (e.g., chr1:1000-2000)</p>
            </div>
            <div className="feature-card">
              <span className="feature-icon">👥</span>
              <h3>Cohort Analysis</h3>
              <p>View all samples in the cohort for the selected region</p>
            </div>
            <div className="feature-card">
              <span className="feature-icon">🔬</span>
              <h3>Sample Details</h3>
              <p>Click on any sample to view detailed information</p>
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="cohort-mode-view">
      <div className="cohort-mode-header">
        <h2>👥 Cohort Analysis</h2>
        <p className="cohort-mode-subtitle">
          Analyze tandem repeat regions across all samples in your cohort
        </p>
      </div>

      {/* Cohort Analysis */}
      {selectedRegion ? (
        <>
          <CohortAnalysis
            mode={mode}
            region={selectedRegion}
            publicVcfFolder={cohortFolderPath}
          />
          
          {/* Floating Navigation for Previous/Next Region */}
          {availableRegions.length > 1 && (
            <div className="floating-navigation">
              <button
                className="floating-nav-btn floating-nav-prev"
                onClick={handlePreviousRegion}
                disabled={!hasPreviousRegion || loadingRegions}
                title="Previous region"
                aria-label="Previous region"
              >
                <span className="floating-nav-icon">◀</span>
              </button>
              
              <div className="floating-nav-info">
                <span className="floating-nav-current">
                  {currentRegionIndex >= 0 ? currentRegionIndex + 1 : '?'}
                </span>
                <span className="floating-nav-separator">/</span>
                <span className="floating-nav-total">{availableRegions.length}</span>
              </div>
              
              <button
                className="floating-nav-btn floating-nav-next"
                onClick={handleNextRegion}
                disabled={!hasNextRegion || loadingRegions}
                title="Next region"
                aria-label="Next region"
              >
                <span className="floating-nav-icon">▶</span>
              </button>
            </div>
          )}
        </>
      ) : (
        <div className="no-region-selected">
          <div className="empty-state-icon">🗺️</div>
          <h3>{cohortFolderPath ? 'No region selected' : 'Cohort folder not loaded'}</h3>
          {cohortFolderPath ? (
            <>
              <p>Please enter a region in the format: <code>chr:start-end</code> (e.g., chr1:1000-2000)</p>
              <p className="hint-text">The region will be analyzed across all samples in your cohort folder.</p>
            </>
          ) : (
            <p>Please load a cohort folder from the sidebar to begin analysis.</p>
          )}
        </div>
      )}
    </div>
  );
};

export default CohortModeView;

