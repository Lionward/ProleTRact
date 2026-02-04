import React, { useEffect, useState, useRef, useMemo } from 'react';
import axios from 'axios';
import './RegionVisualization.css';
import PopulationComparison from './PopulationComparison';
import PathogenicityPanel from './PathogenicityPanel';
import AnnotationPanel from './AnnotationPanel';
import RegionInfoCard from './RegionInfoCard';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

interface Record {
  chr: string;
  pos: number;
  stop: number;
  motifs: string[];
  motif_ids_h1: string[];
  motif_ids_h2: string[];
  motif_ids_ref: string[];
  ref_CN: number | null;
  CN_H1: string | null;
  CN_H2: string | null;
  spans: string[];
  ref_allele: string;
  alt_allele1: string;
  alt_allele2: string;
  gt: string;
  supported_reads_h1: number;
  supported_reads_h2: number;
  id: string;
}

interface RegionVisualizationProps {
  region: string;
  vcfPath: string;
  publicVcfFolder?: string;
  onRegionSelect?: (region: string) => void;
}

const RegionVisualization: React.FC<RegionVisualizationProps> = ({ 
  region, 
  vcfPath,
  publicVcfFolder,
  onRegionSelect
}) => {
  const [record, setRecord] = useState<Record | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');
  const [availableRegions, setAvailableRegions] = useState<string[]>([]);
  const [activeViewTab, setActiveViewTab] = useState<'sequence' | 'population'>('sequence');
  const [pathogenicInfo, setPathogenicInfo] = useState<{ 
    pathogenic_threshold?: number;
    gene?: string;
    disease?: string;
    inheritance?: string;
  } | null>(null);
  const [viewMode, setViewMode] = useState<ViewMode>('sequence');
  const scrollPositionRef = useRef<number>(0);

  // get regions for the search dropdown
  useEffect(() => {
    if (!vcfPath) {
      setAvailableRegions([]);
      return;
    }

    const fetchRegions = async () => {
      try {
        const response = await axios.get(`${API_BASE}/api/vcf/regions`, {
          params: { vcf_path: vcfPath }
        });
        setAvailableRegions(response.data.regions || []);
      } catch (err) {
        // if this fails, autocomplete just wont work, no big deal
        console.error('Failed to fetch regions for autocomplete:', err);
      }
    };

    fetchRegions();
  }, [vcfPath]);

  // keep track of scroll position as user scrolls
  useEffect(() => {
    const handleScroll = () => {
      scrollPositionRef.current = window.scrollY || document.documentElement.scrollTop;
    };

    window.addEventListener('scroll', handleScroll, { passive: true });
    return () => window.removeEventListener('scroll', handleScroll);
  }, []);

  useEffect(() => {
    if (!region || !vcfPath) {
      setRecord(null);
      return;
    }

    // save where we are scrolled to before switching regions
    const savedScrollPosition = scrollPositionRef.current || (window.scrollY || document.documentElement.scrollTop);

    const fetchRegionData = async () => {
      setLoading(true);
      setError(null);
      try {
        const response = await axios.get(`${API_BASE}/api/vcf/region/${encodeURIComponent(region)}`, {
          params: { vcf_path: vcfPath }
        });
        setRecord(response.data.record);
      } catch (err: any) {
        setError(err.response?.data?.detail || 'Failed to load region data');
        setRecord(null);
      } finally {
        setLoading(false);
        // restore scroll position - try a few times since the DOM might not be ready yet
        const restoreScroll = () => {
          window.scrollTo({
            top: savedScrollPosition,
            behavior: 'auto'
          });
        };
        // try right away
        restoreScroll();
        // try again after the next frame
        requestAnimationFrame(restoreScroll);
        // and once more after a tiny delay just in case
        setTimeout(restoreScroll, 10);
      }
    };

    fetchRegionData();
  }, [region, vcfPath]);

  // Fetch pathogenic info for the region
  useEffect(() => {
    const fetchPathogenicInfo = async () => {
      if (!region || !record) {
        setPathogenicInfo(null);
        return;
      }
      
      // Parse region to get chr, start, end
      const regionMatch = region.match(/^([^:]+):(\d+)-(\d+)$/);
      if (!regionMatch) {
        setPathogenicInfo(null);
        return;
      }
      
      const [, chr, startStr, endStr] = regionMatch;
      const start = parseInt(startStr, 10);
      const end = parseInt(endStr, 10);
      
      try {
        const response = await axios.get(`${API_BASE}/api/pathogenic/check`, {
          params: {
            chr: chr,
            start: start,
            end: end
          }
        });
        
        if (response.data.pathogenic) {
          const threshold = response.data.pathogenic_threshold;
          setPathogenicInfo({ 
            pathogenic_threshold: threshold !== undefined && threshold !== null ? threshold : undefined,
            gene: response.data.gene || undefined,
            disease: response.data.disease || undefined,
            inheritance: response.data.inheritance || undefined
          });
        } else {
          setPathogenicInfo(null);
        }
      } catch (err: any) {
        // Silently fail if endpoint doesn't exist yet
        if (err.response?.status !== 404) {
          console.error('Error checking pathogenicity:', err);
        }
        setPathogenicInfo(null);
      }
    };
    
    fetchPathogenicInfo();
  }, [region, record]);

  const handleSearch = (query: string) => {
    setSearchQuery(query);
  };

  const handleSearchSubmit = async (selectedRegion?: string) => {
    const regionToSearch = selectedRegion || searchQuery;
    if (!regionToSearch || !onRegionSelect) return;
    
    const trimmedQuery = regionToSearch.trim();
    
    // Try to parse and validate the region format
    const regionPattern = /^([\w]+):(\d+)-(\d+)$/;
    const match = trimmedQuery.match(regionPattern);
    
    if (match) {
      // It's a valid region format
      onRegionSelect(trimmedQuery);
      setSearchQuery('');
    } else {
      // It might be a gene name - search for it
      try {
        const response = await axios.get(`${API_BASE}/api/pathogenic/search`, {
          params: { gene: trimmedQuery }
        });
        if (response.data.success && response.data.regions && response.data.regions.length > 0) {
          const firstRegion = response.data.regions[0];
          onRegionSelect(firstRegion.region);
          setSearchQuery('');
        } else {
          // Show error or message that no gene found
          console.log(`No regions found for gene: ${trimmedQuery}`);
        }
      } catch (err) {
        console.error('Error searching by gene:', err);
      }
    }
  };

  // Extract current alleles from record
  const currentAlleles = useMemo(() => {
    if (!record) return [];
    const alleles: string[] = [];
    if (record.alt_allele1) alleles.push(record.alt_allele1);
    if (record.alt_allele2) alleles.push(record.alt_allele2);
    return alleles;
  }, [record]);

  if (loading) {
    return <div className="loading-message">Loading region data...</div>;
  }

  if (error) {
    return <div className="error-message">Error: {error}</div>;
  }

  if (!record) {
    return (
      <div className="region-visualization-container">
        <div style={{ marginBottom: '1rem' }}>
          <RegionSearchBar 
            value={searchQuery} 
            onChange={handleSearch}
            onSubmit={handleSearchSubmit}
            placeholder="Search region (e.g., chr1:1000-2000) or gene name"
            availableRegions={availableRegions}
          />
        </div>
        <div className="no-data-message">No data available for this region</div>
      </div>
    );
  }

  return (
    <div className="region-visualization-container">
      <div style={{ marginBottom: '1rem' }}>
        <RegionSearchBar 
          value={searchQuery} 
          onChange={handleSearch}
          onSubmit={handleSearchSubmit}
          placeholder={`Search region (e.g., ${record.chr}:${record.pos}-${record.stop}) or gene name`}
          availableRegions={availableRegions}
        />
      </div>
      
      {/* Tabs for Sequence View and Population Comparison */}
      <div className="region-view-tabs">
        <button
          className={`view-tab-button ${activeViewTab === 'sequence' ? 'active' : ''}`}
          onClick={() => setActiveViewTab('sequence')}
        >
          🔬 Sequence View
        </button>
        <button
          className={`view-tab-button ${activeViewTab === 'population' ? 'active' : ''}`}
          onClick={() => setActiveViewTab('population')}
        >
          🌐 Allele vs Population
        </button>
      </div>

      {/* Keep all tabs mounted but hidden - prevents recalculation on tab switch */}
      <div className={`view-tab-content ${activeViewTab === 'sequence' ? 'tab-active' : 'tab-hidden'}`}>
        <RegionInfoCard
          chr={record.chr}
          pos={record.pos}
          stop={record.stop}
          region={region}
        />
        <PathogenicityPanel
          chr={record.chr}
          pos={record.pos}
          stop={record.stop}
          region={region}
          record={record}
        />
        <AnnotationPanel region={region} />
        <RegionDisplay 
          record={record} 
          pathogenicThreshold={pathogenicInfo?.pathogenic_threshold}
          viewMode={viewMode}
          onViewModeChange={setViewMode}
        />
      </div>

      <div className={`view-tab-content ${activeViewTab === 'population' ? 'tab-active' : 'tab-hidden'}`}>
        <RegionInfoCard
          chr={record.chr}
          pos={record.pos}
          stop={record.stop}
          region={region}
        />
        <PopulationComparison 
          record={record} 
          region={region}
          publicVcfFolder={publicVcfFolder}
        />
      </div>
    </div>
  );
};

interface RegionSearchBarProps {
  value: string;
  onChange: (value: string) => void;
  onSubmit: (selectedRegion?: string) => void;
  placeholder: string;
  availableRegions: string[];
  showAllWhenEmpty?: boolean; // Show all regions when input is empty (for filter mode)
}

const RegionSearchBar: React.FC<RegionSearchBarProps> = ({ 
  value, 
  onChange, 
  onSubmit, 
  placeholder,
  availableRegions,
  showAllWhenEmpty = false
}) => {
  const [showSuggestions, setShowSuggestions] = useState(false);
  const [highlightedIndex, setHighlightedIndex] = useState(-1);
  const [geneSearchResults, setGeneSearchResults] = useState<Array<{
    region: string;
    chr: string;
    start: number;
    end: number;
    gene: string;
    disease?: string;
    inheritance?: string;
    motif?: string;
    pathogenic_threshold?: number;
  }>>([]);
  const [isSearchingGene, setIsSearchingGene] = useState(false);
  const searchRef = useRef<HTMLDivElement>(null);
  const inputRef = useRef<HTMLInputElement>(null);

  // Check if query is a region format (chr:start-end) or a gene name
  const isRegionFormat = (query: string): boolean => {
    const regionPattern = /^([\w]+):(\d+)-(\d+)$/;
    return regionPattern.test(query.trim());
  };

  // Filter regions based on search query
  const filteredRegions = React.useMemo(() => {
    const query = value ? value.toLowerCase().trim() : '';
    
    // If no query and showAllWhenEmpty is true, show all available regions
    if (!query && showAllWhenEmpty) {
      // Show all available regions (which are already filtered if pathogenic filter is on)
      return availableRegions.slice(0, 50); // Show up to 50 when filter is active
    }
    
    // If no query and showAllWhenEmpty is false, return empty (original behavior)
    if (!query) {
      return [];
    }
    
    // If it's a region format, don't filter by substring
    if (isRegionFormat(value)) {
      return availableRegions.filter(region => region === value.trim()).slice(0, 1);
    }
    
    return availableRegions
      .filter(region => region.toLowerCase().includes(query))
      .slice(0, 10); // Limit to 10 suggestions
  }, [value, availableRegions, showAllWhenEmpty]);

  // Search by gene when query is not a region format
  useEffect(() => {
    const query = value.trim();
    if (!query || isRegionFormat(query)) {
      setGeneSearchResults([]);
      setIsSearchingGene(false);
      return;
    }

    // Debounce gene search
    const timeoutId = setTimeout(async () => {
      setIsSearchingGene(true);
      try {
        const response = await axios.get(`${API_BASE}/api/pathogenic/search`, {
          params: { gene: query }
        });
        if (response.data.success && response.data.regions) {
          setGeneSearchResults(response.data.regions.slice(0, 10)); // Limit to 10 results
        } else {
          setGeneSearchResults([]);
        }
      } catch (err) {
        console.error('Error searching by gene:', err);
        setGeneSearchResults([]);
      } finally {
        setIsSearchingGene(false);
      }
    }, 300); // 300ms debounce

    return () => clearTimeout(timeoutId);
  }, [value]);

  // Handle click outside to close suggestions
  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (searchRef.current && !searchRef.current.contains(event.target as Node)) {
        setShowSuggestions(false);
      }
    };

    document.addEventListener('mousedown', handleClickOutside);
    return () => {
      document.removeEventListener('mousedown', handleClickOutside);
    };
  }, []);

  const handleInputChange = (newValue: string) => {
    onChange(newValue);
    setShowSuggestions(true);
    setHighlightedIndex(-1);
  };

  const handleSelectRegion = (selectedRegion: string) => {
    onChange(selectedRegion);
    setShowSuggestions(false);
    onSubmit(selectedRegion);
  };

  const allSuggestions = React.useMemo(() => {
    const suggestions: Array<{ type: 'region' | 'gene'; value: string; label: string; data?: any }> = [];
    
    // Add region suggestions
    filteredRegions.forEach(region => {
      suggestions.push({ type: 'region', value: region, label: region });
    });
    
    // Add gene search results
    geneSearchResults.forEach(result => {
      const label = `${result.region}${result.gene ? ` (${result.gene})` : ''}${result.disease ? ` - ${result.disease}` : ''}`;
      suggestions.push({ 
        type: 'gene', 
        value: result.region, 
        label,
        data: result
      });
    });
    
    return suggestions;
  }, [filteredRegions, geneSearchResults]);

  // Auto-show suggestions when filter is active and there are filtered regions
  useEffect(() => {
    if (showAllWhenEmpty && filteredRegions.length > 0) {
      setShowSuggestions(true);
    }
  }, [showAllWhenEmpty, filteredRegions.length]);

  const handleKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (!showSuggestions || allSuggestions.length === 0) {
      if (e.key === 'Enter') {
        onSubmit();
      }
      return;
    }

    switch (e.key) {
      case 'ArrowDown':
        e.preventDefault();
        setHighlightedIndex(prev => 
          prev < allSuggestions.length - 1 ? prev + 1 : prev
        );
        break;
      case 'ArrowUp':
        e.preventDefault();
        setHighlightedIndex(prev => prev > 0 ? prev - 1 : -1);
        break;
      case 'Enter':
        e.preventDefault();
        if (highlightedIndex >= 0 && highlightedIndex < allSuggestions.length) {
          handleSelectRegion(allSuggestions[highlightedIndex].value);
        } else {
          onSubmit();
        }
        break;
      case 'Escape':
        setShowSuggestions(false);
        setHighlightedIndex(-1);
        break;
    }
  };

  return (
    <div className="region-search-bar-container" ref={searchRef}>
      <div className="region-search-bar">
        <div className="search-icon">🔍</div>
        <input
          ref={inputRef}
          type="text"
          className="search-input"
          value={value}
          onChange={(e) => handleInputChange(e.target.value)}
          onFocus={() => setShowSuggestions(true)}
          onKeyDown={handleKeyDown}
          placeholder={placeholder}
        />
        {value && (
          <button 
            className="search-button"
            onClick={() => onSubmit()}
            title="Search region"
          >
            Go
          </button>
        )}
      </div>
      {showSuggestions && (allSuggestions.length > 0 || isSearchingGene) && (
        <div className="search-suggestions">
          {isSearchingGene && allSuggestions.length === 0 && (
            <div className="suggestion-item" style={{ fontStyle: 'italic', color: '#9ca3af' }}>
              Searching for gene...
            </div>
          )}
          {allSuggestions.map((suggestion, index) => (
            <div
              key={`${suggestion.type}-${suggestion.value}-${index}`}
              className={`suggestion-item ${index === highlightedIndex ? 'highlighted' : ''}`}
              onClick={() => handleSelectRegion(suggestion.value)}
              onMouseEnter={() => setHighlightedIndex(index)}
              title={suggestion.type === 'gene' && suggestion.data?.disease ? suggestion.data.disease : undefined}
            >
              {suggestion.type === 'gene' && <span style={{ marginRight: '8px' }}>🧬</span>}
              {suggestion.label}
            </div>
          ))}
          {!isSearchingGene && allSuggestions.length === 0 && value.trim() && (
            <div className="suggestion-item" style={{ fontStyle: 'italic', color: '#9ca3af' }}>
              No matches found. Try a region (e.g., chr1:1000-2000) or gene name.
            </div>
          )}
        </div>
      )}
    </div>
  );
};

// Color palette matching ProleTRact
const getColorPalette = (numColors: number): string[] => {
  const primaryPalette = [
    '#667eea', '#764ba2', '#4fd1c7', '#68d391',
    '#f6ad55', '#fc8181', '#f093fb', '#f6e05e'
  ];
  const secondaryPalette = [
    '#5a67d8', '#6b46c1', '#38a169', '#dd6b20',
    '#e53e3e', '#d53f8c', '#d69e2e', '#319795'
  ];
  const combined = [...primaryPalette, ...secondaryPalette];
  
  const goldenRatio = 0.618033988749895;
  const colors = [...combined];
  for (let i = combined.length; i < numColors; i++) {
    const hue = (i * goldenRatio) % 1;
    const saturation = 0.6 + 0.2 * (i % 3) / 3;
    const value = 0.7 + 0.2 * ((i + 1) % 2);
    const rgb = hsvToRgb(hue, saturation, value);
    const hex = `#${rgb[0].toString(16).padStart(2, '0')}${rgb[1].toString(16).padStart(2, '0')}${rgb[2].toString(16).padStart(2, '0')}`;
    colors.push(hex);
  }
  
  return colors.slice(0, numColors);
};

const hsvToRgb = (h: number, s: number, v: number): [number, number, number] => {
  const i = Math.floor(h * 6);
  const f = h * 6 - i;
  const p = v * (1 - s);
  const q = v * (1 - f * s);
  const t = v * (1 - (1 - f) * s);
  
  switch (i % 6) {
    case 0: return [Math.round(v * 255), Math.round(t * 255), Math.round(p * 255)];
    case 1: return [Math.round(q * 255), Math.round(v * 255), Math.round(p * 255)];
    case 2: return [Math.round(p * 255), Math.round(v * 255), Math.round(t * 255)];
    case 3: return [Math.round(p * 255), Math.round(q * 255), Math.round(v * 255)];
    case 4: return [Math.round(t * 255), Math.round(p * 255), Math.round(v * 255)];
    default: return [Math.round(v * 255), Math.round(p * 255), Math.round(q * 255)];
  }
};

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

const calculateMotifCoverage = (sequence: string, spans: string): number => {
  if (!sequence || !spans) return 0;
  const ranges = parseMotifRange(spans);
  let motifLength = 0;
  ranges.forEach(([start, end]) => {
    motifLength += (end - start + 1);
  });
  return sequence.length > 0 ? Math.round((motifLength / sequence.length) * 100) : 0;
};

const calculateMotifCount = (motifIds: string[]): number => {
  return motifIds.filter(id => id !== '.' && id !== '').length;
};

interface SequenceMetadata {
  length: number;
  motifCoverage: number;
  motifCount: number;
  supportingReads?: number;
}

type ViewMode = 'sequence' | 'bar';

interface SequenceTrackProps {
  sequence: string;
  motifIds: string[];
  spans: string;
  label: string;
  metadata: SequenceMetadata;
  motifColors: { [key: number]: string };
  motifNames: string[];
  isReference?: boolean;
  isPathogenic?: boolean;
  hasPathogenicThreshold?: boolean; // Whether a pathogenic threshold is defined for this region
  viewMode?: ViewMode; // 'sequence' or 'bar'
}

const SequenceTrack: React.FC<SequenceTrackProps> = ({
  sequence,
  motifIds,
  spans,
  label,
  metadata,
  motifColors,
  motifNames,
  isReference = false,
  isPathogenic = false,
  hasPathogenicThreshold = false,
  viewMode = 'sequence'
}) => {
  const scrollContainerRef = useRef<HTMLDivElement>(null);
  const sequenceRef = useRef<HTMLDivElement>(null);
  
  if (!sequence || sequence === '.' || sequence === '') {
    return (
      <div className={`sequence-track ${isReference ? 'reference-track' : ''}`}>
        <div className="sequence-header">
          <div className="sequence-label-container">
            <span className="sequence-label">{label}</span>
          </div>
          <SequenceMetadataDisplay metadata={metadata} />
        </div>
        <div className="sequence-content no-sequence">No sequence data</div>
      </div>
    );
  }

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

  // Calculate scale markers (every 25% and endpoints)
  const sequenceLength = sequence.length;
  const scaleMarkers = [0];
  if (sequenceLength > 0) {
    scaleMarkers.push(Math.round(sequenceLength * 0.25));
    scaleMarkers.push(Math.round(sequenceLength * 0.5));
    scaleMarkers.push(Math.round(sequenceLength * 0.75));
    scaleMarkers.push(sequenceLength);
  }
  const uniqueMarkers = Array.from(new Set(scaleMarkers)).sort((a, b) => a - b);

  // Character width approximation (monospace font)
  const charWidth = 8.4; // Approximate width for monospace font at 0.875rem
  const sequenceWidth = sequenceLength * charWidth;

  // Only apply color coding if a pathogenic threshold is defined
  const shouldApplyStyling = hasPathogenicThreshold && !isReference;
  
  // Render bar view
  if (viewMode === 'bar') {
    
    return (
      <div 
        className={`sequence-track bar-view ${isReference ? 'reference-track' : ''} ${shouldApplyStyling ? (isPathogenic ? 'pathogenic-allele' : 'non-pathogenic-allele') : ''}`}
        style={shouldApplyStyling ? (isPathogenic ? {
          border: '3px solid #F87171',
          borderRadius: '6px',
          backgroundColor: 'rgba(248, 113, 113, 0.05)',
          boxShadow: '0 0 8px rgba(248, 113, 113, 0.3)'
        } : {
          border: '3px solid #10B981',
          borderRadius: '6px',
          backgroundColor: 'rgba(16, 185, 129, 0.05)',
          boxShadow: '0 0 8px rgba(16, 185, 129, 0.2)'
        }) : {}}
      >
        <div className="sequence-header">
          <div className="sequence-label-container">
            <span className="sequence-label">{label}</span>
          </div>
          <SequenceMetadataDisplay metadata={metadata} />
        </div>
        
        <div className="bar-view-container">
          <div className="bar-view-wrapper">
            {segments.map((seg, idx) => {
              const segLength = seg.end - seg.start + 1;
              const segWidthPercent = (segLength / sequenceLength) * 100;
              
              if (seg.type === 'interruption') {
                return (
                  <div
                    key={idx}
                    className="bar-segment interruption"
                    style={{ width: `${segWidthPercent}%` }}
                    title={`Interruption\nPosition: ${seg.start + 1}-${seg.end + 1}\nLength: ${segLength} bp`}
                  />
                );
              } else {
                const color = seg.motifId !== undefined ? motifColors[seg.motifId] || '#cccccc' : '#cccccc';
                const motifName = seg.motifId !== undefined && motifNames[seg.motifId] 
                  ? motifNames[seg.motifId] 
                  : `Motif ${seg.motifId}`;
                return (
                  <div
                    key={idx}
                    className="bar-segment motif"
                    style={{ 
                      width: `${segWidthPercent}%`,
                      backgroundColor: color,
                      borderColor: `rgba(0, 0, 0, 0.4)`
                    }}
                    title={`${motifName}\nPosition: ${seg.start + 1}-${seg.end + 1}\nLength: ${segLength} bp`}
                  />
                );
              }
            })}
          </div>
          
          <div className="bar-view-scale">
            <div className="scale-line" />
            {uniqueMarkers.map((pos, idx) => (
              <div 
                key={idx}
                className="scale-marker"
                style={{ left: `${(pos / sequenceLength) * 100}%` }}
              >
                <div className="scale-tick" />
                <div className="scale-label">{pos}</div>
              </div>
            ))}
          </div>
        </div>
      </div>
    );
  }
  
  // Render sequence view (default)
  return (
    <div 
      className={`sequence-track ${isReference ? 'reference-track' : ''} ${shouldApplyStyling ? (isPathogenic ? 'pathogenic-allele' : 'non-pathogenic-allele') : ''}`}
      style={shouldApplyStyling ? (isPathogenic ? {
        border: '3px solid #F87171', // Softer red/coral color
        borderRadius: '6px',
        backgroundColor: 'rgba(248, 113, 113, 0.05)', // Very light red background
        boxShadow: '0 0 8px rgba(248, 113, 113, 0.3)'
      } : {
        border: '3px solid #10B981', // Green color
        borderRadius: '6px',
        backgroundColor: 'rgba(16, 185, 129, 0.05)', // Very light green background
        boxShadow: '0 0 8px rgba(16, 185, 129, 0.2)'
      }) : {}}
    >
      <div className="sequence-header">
        <div className="sequence-label-container">
          <span className="sequence-label">{label}</span>
        </div>
        <SequenceMetadataDisplay metadata={metadata} />
      </div>
      
      <div 
        className="sequence-scroll-container"
        ref={scrollContainerRef}
      >
        <div 
          className="sequence-content-wrapper"
          style={{ minWidth: `${Math.max(sequenceWidth, 100)}px` }}
        >
          <div className="sequence-line" ref={sequenceRef}>
            {segments.map((seg, idx) => {
              const seq = sequence.substring(seg.start, seg.end + 1);
              if (seg.type === 'interruption') {
                return (
                  <span 
                    key={idx} 
                    className="sequence-segment interruption" 
                    title={`Interruption: ${seq}\nPosition: ${seg.start + 1}-${seg.end + 1}`}
                  >
                    {seq}
                  </span>
                );
              } else {
                const color = seg.motifId !== undefined ? motifColors[seg.motifId] || '#cccccc' : '#cccccc';
                const motifName = seg.motifId !== undefined && motifNames[seg.motifId] 
                  ? motifNames[seg.motifId] 
                  : `Motif ${seg.motifId}`;
                return (
                  <span
                    key={idx}
                    className="sequence-segment motif"
                    style={{ 
                      backgroundColor: color,
                      color: '#1f2937',
                      borderColor: 'rgba(0, 0, 0, 0.25)'
                    }}
                    title={`${motifName}\nSequence: ${seq}\nPosition: ${seg.start + 1}-${seg.end + 1}\nLength: ${seg.end - seg.start + 1} bp`}
                  >
                    {seq}
                  </span>
                );
              }
            })}
          </div>
          
          <div className="sequence-scale">
            <div className="scale-line" />
            {uniqueMarkers.map((pos, idx) => (
              <div 
                key={idx}
                className="scale-marker"
                style={{ left: `${(pos / sequenceLength) * 100}%` }}
              >
                <div className="scale-tick" />
                <div className="scale-label">{pos}</div>
              </div>
            ))}
          </div>
        </div>
      </div>
    </div>
  );
};

interface SequenceMetadataDisplayProps {
  metadata: SequenceMetadata;
}

const SequenceMetadataDisplay: React.FC<SequenceMetadataDisplayProps> = ({ metadata }) => {
  return (
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
  );
};

interface RegionDisplayProps {
  record: Record;
  pathogenicThreshold?: number;
  viewMode: ViewMode;
  onViewModeChange: (mode: ViewMode) => void;
}

const RegionDisplay: React.FC<RegionDisplayProps> = ({ record, pathogenicThreshold, viewMode, onViewModeChange }) => {
  
  const motifColors = getColorPalette(record.motifs.length);
  const motifColorMap: { [key: number]: string } = {};
  record.motifs.forEach((_, idx) => {
    motifColorMap[idx] = motifColors[idx];
  });

  const refMetadata: SequenceMetadata = {
    length: record.ref_allele?.length || 0,
    motifCoverage: calculateMotifCoverage(record.ref_allele, record.spans[0] || ''),
    motifCount: calculateMotifCount(record.motif_ids_ref),
  };

  const h1Metadata: SequenceMetadata = {
    length: record.alt_allele1?.length || 0,
    motifCoverage: calculateMotifCoverage(record.alt_allele1, record.spans[1] || ''),
    motifCount: calculateMotifCount(record.motif_ids_h1),
    supportingReads: record.supported_reads_h1,
  };

  const h2Metadata: SequenceMetadata = {
    length: record.alt_allele2?.length || 0,
    motifCoverage: calculateMotifCoverage(record.alt_allele2, record.spans[2] || ''),
    motifCount: calculateMotifCount(record.motif_ids_h2),
    supportingReads: record.supported_reads_h2,
  };

  return (
    <div className="region-visualization">
      <div className="region-header">
        <div className="region-title-section">
          <h2 className="region-title">Region: {record.chr}:{record.pos}-{record.stop}</h2>
          <div className="region-id">ID: {record.id}</div>
        </div>
        <div className="region-stats">
          <div className="stat-card">
            <div className="stat-label">Genotype</div>
            <div className="stat-value">{record.gt}</div>
          </div>
          <div className="stat-card">
            <div className="stat-label">Ref CN</div>
            <div className="stat-value">{record.ref_CN ?? '—'}</div>
          </div>
          <div className="stat-card">
            <div className="stat-label">H1 CN</div>
            <div className="stat-value">{record.CN_H1 ?? '—'}</div>
          </div>
          <div className="stat-card">
            <div className="stat-label">H2 CN</div>
            <div className="stat-value">{record.CN_H2 ?? '—'}</div>
          </div>
        </div>
      </div>

      <div className="motif-legend">
        <div className="legend-header">Motif Legend</div>
        <div className="legend-items">
          {record.motifs.map((motif, idx) => (
            <div key={idx} className="legend-item">
              <span
                className="legend-color"
                style={{ backgroundColor: motifColorMap[idx] }}
              />
              <span className="legend-text">{motif} ({motif.length}bp)</span>
            </div>
          ))}
          <div className="legend-item">
            <span className="legend-color interruption" />
            <span className="legend-text">Interruption</span>
          </div>
        </div>
      </div>

      {/* View Mode Controls */}
      <div className="view-mode-controls">
        <div className="view-mode-group">
          <label className="view-mode-label">Display Mode:</label>
          <div className="view-mode-buttons">
            <button
              className={`view-mode-btn ${viewMode === 'sequence' ? 'active' : ''}`}
              onClick={() => onViewModeChange('sequence')}
              title="Show actual DNA sequence with color-coded motifs"
            >
              📝 Sequence
            </button>
            <button
              className={`view-mode-btn ${viewMode === 'bar' ? 'active' : ''}`}
              onClick={() => onViewModeChange('bar')}
              title="Show motifs as colored bars (absolute scale - shows actual differences between haplotypes)"
            >
              📊 Bar
            </button>
          </div>
        </div>
      </div>

      <div className="sequences-container">
        <SequenceTrack
          sequence={record.ref_allele}
          motifIds={record.motif_ids_ref}
          spans={record.spans[0] || ''}
          label="Reference"
          metadata={refMetadata}
          motifColors={motifColorMap}
          motifNames={record.motifs}
          isReference={true}
          viewMode={viewMode}
        />
        {/* Show Allele 1 - either as sequence track or deleted indicator */}
        {(() => {
          const cnH1 = record.CN_H1;
          const isDeleted = cnH1 === '0' || cnH1 === null || cnH1 === '';
          const hasSequence = record.alt_allele1 && record.alt_allele1 !== '' && record.alt_allele1 !== '.';
          
          if (isDeleted) {
            return (
              <div className="sequence-track deleted-allele">
                <div className="sequence-header">
                  <div className="sequence-label-container">
                    <span className="sequence-label">Allele 1</span>
                    <span className="deleted-badge">Deleted</span>
                  </div>
                  <div className="sequence-metadata">
                    <div className="metadata-item">
                      <span className="metadata-label">Copy Number:</span>
                      <span className="metadata-value">0</span>
                    </div>
                  </div>
                </div>
                <div className="sequence-content deleted-content">
                  <div className="deleted-message">This allele has been deleted (CN = 0)</div>
                </div>
              </div>
            );
          } else if (hasSequence) {
            // Check if this allele exceeds pathogenic threshold
            const hasThreshold = pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0;
            const motifCount = calculateMotifCount(record.motif_ids_h1);
            const isPathogenic = hasThreshold && motifCount >= pathogenicThreshold;
            
            return (
              <SequenceTrack
                sequence={record.alt_allele1}
                motifIds={record.motif_ids_h1}
                spans={record.spans[1] || ''}
                label="Allele 1"
                metadata={h1Metadata}
                motifColors={motifColorMap}
                motifNames={record.motifs}
                isPathogenic={isPathogenic}
                hasPathogenicThreshold={hasThreshold}
                viewMode={viewMode}
              />
            );
          }
          return null;
        })()}
        {/* Show Allele 2 - either as sequence track or deleted indicator */}
        {(() => {
          const cnH2 = record.CN_H2;
          const isDeleted = cnH2 === '0' || cnH2 === null || cnH2 === '';
          const hasSequence = record.alt_allele2 && record.alt_allele2 !== '' && record.alt_allele2 !== '.';
          
          if (isDeleted) {
            return (
              <div className="sequence-track deleted-allele">
                <div className="sequence-header">
                  <div className="sequence-label-container">
                    <span className="sequence-label">Allele 2</span>
                    <span className="deleted-badge">Deleted</span>
                  </div>
                  <div className="sequence-metadata">
                    <div className="metadata-item">
                      <span className="metadata-label">Copy Number:</span>
                      <span className="metadata-value">0</span>
                    </div>
                  </div>
                </div>
                <div className="sequence-content deleted-content">
                  <div className="deleted-message">This allele has been deleted (CN = 0)</div>
                </div>
              </div>
            );
          } else if (hasSequence) {
            // Check if this allele exceeds pathogenic threshold
            const hasThreshold = pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0;
            const motifCount = calculateMotifCount(record.motif_ids_h2);
            const isPathogenic = hasThreshold && motifCount >= pathogenicThreshold;
            
            return (
              <SequenceTrack
                sequence={record.alt_allele2}
                motifIds={record.motif_ids_h2}
                spans={record.spans[2] || ''}
                label="Allele 2"
                metadata={h2Metadata}
                motifColors={motifColorMap}
                motifNames={record.motifs}
                isPathogenic={isPathogenic}
                hasPathogenicThreshold={hasThreshold}
                viewMode={viewMode}
              />
            );
          }
          return null;
        })()}
      </div>
    </div>
  );
};

export default RegionVisualization;
export { RegionSearchBar };
