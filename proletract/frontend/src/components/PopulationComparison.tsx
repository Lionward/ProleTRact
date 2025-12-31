import React, { useState, useEffect, useMemo, memo, useRef, useCallback } from 'react';
import axios from 'axios';
import './PopulationComparison.css';
import ExportMenu from './ExportMenu';
import PopulationFrequencyPanel from './PopulationFrequencyPanel';
import { exportToFASTA, generateFilename } from '../utils/exportUtils';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

// Color palette for motifs (defined early for use in helper functions)
const COLOR_PALETTE = [
  '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
  '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6',
  '#60a5fa', '#34d399', '#fbbf24', '#fb7185', '#a78bfa', '#38bdf8'
];

// Cache for processed data by region
const processedDataCache = new Map<string, {
  genotypeComparisonData: { [sample: string]: string };
  sequencesData: { sequences: Array<{ name: string; sequence: string }>; spanList: string[]; motifIdsList: string[][]; sortedSamples: string[] };
  populationRecords: { [key: string]: PopulationRecord };
}>();

// Cache for combined plot data
const combinedPlotCache = new Map<string, {
  motifData: Array<{ sample: string; start: number; end: number; motif: string; motifIndex: number; length: number }>;
  heatmapData: Array<{ sample: string; motif: string; count: number }>;
  motifColorMap: { [key: string]: string };
  maxLength: number;
  maxCount: number;
  sortedMotifs: string[];
  sampleDataMap: { [sample: string]: Array<{ start: number; end: number; motif: string; motifIndex: number; length: number }> };
}>();

// Cache for parsed motif ranges to avoid re-parsing
const parsedRangeCache = new Map<string, Array<[number, number]>>();

interface PopulationRecord {
  chr: string;
  pos: number;
  stop: number;
  motifs: string[];
  motif_ids_h: string[];
  motif_ids_ref: string[];
  ref_CN: number;
  CN_H: number;
  spans: string;
  ref_allele: string;
  alt_allele: string;
  gt: string;
  id: string;
}

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

interface PopulationComparisonProps {
  record: Record;
  region: string;
  publicVcfFolder?: string;
}

const PopulationComparison: React.FC<PopulationComparisonProps> = ({ record, region, publicVcfFolder }) => {
  const [populationRecords, setPopulationRecords] = useState<{ [key: string]: PopulationRecord }>({});
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [pathogenicInfo, setPathogenicInfo] = useState<{ 
    pathogenic_threshold?: number;
    gene?: string;
    disease?: string;
    inheritance?: string;
  } | null>(null);
  
  // Collapsible sections state
  const [expandedSections, setExpandedSections] = useState({
    frequency: true,
    genotype: true,
    stackHeatmap: true,
    barPlot: true,
    cluster: true
  });
  
  // Virtual scrolling state
  const [virtualScrollEnabled, setVirtualScrollEnabled] = useState(false);
  const [visibleRange, setVisibleRange] = useState({ start: 0, end: 50 });
  const VIRTUAL_SCROLL_THRESHOLD = 50; // Enable virtual scrolling if more than 50 samples
  
  // Search and filter state
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedGenotypes, setSelectedGenotypes] = useState<string[]>([]);
  
  // Sort state
  const [sortBy, setSortBy] = useState<'length' | 'motifCount' | 'name' | 'genotype'>('length');
  const [sortDirection, setSortDirection] = useState<'asc' | 'desc'>('desc');

  // Create cache key
  const cacheKey = useMemo(() => {
    if (!publicVcfFolder || !region) return null;
    return `${publicVcfFolder}:${region}:${record.id}`;
  }, [publicVcfFolder, region, record.id]);

  // Try to get cached data immediately - this makes switching instant
  const cachedData = useMemo(() => {
    if (!cacheKey) return null;
    return processedDataCache.get(cacheKey) || null;
  }, [cacheKey]);

  // Use cached data if available, otherwise use state
  const effectivePopulationRecords = useMemo(() => {
    if (cachedData && Object.keys(cachedData.populationRecords).length > 0) {
      return cachedData.populationRecords;
    }
    return populationRecords;
  }, [cachedData, populationRecords]);

  // Process population records and create genotype comparison data
  // Must be called before any early returns (React Hooks rules)
  const genotypeComparisonData = useMemo(() => {
    if (!publicVcfFolder) return {};
    
    // Check cache first - instant return
    if (cachedData && Object.keys(cachedData.genotypeComparisonData).length > 0) {
      return cachedData.genotypeComparisonData;
    }
    
    if (Object.keys(effectivePopulationRecords).length === 0) {
      return {};
    }
    
    const data = computeGenotypeComparison(record, effectivePopulationRecords);
    
    // Cache the result
    if (cacheKey) {
      const cached = processedDataCache.get(cacheKey);
      if (cached) {
        cached.genotypeComparisonData = data;
      } else {
        processedDataCache.set(cacheKey, {
          genotypeComparisonData: data,
          sequencesData: { sequences: [], spanList: [], motifIdsList: [], sortedSamples: [] },
          populationRecords: effectivePopulationRecords
        });
      }
    }
    
    return data;
  }, [record, effectivePopulationRecords, publicVcfFolder, cacheKey, cachedData]);

  // Create sequences data for stack plot and bar plot
  // Must be called before any early returns (React Hooks rules)
  const sequencesData = useMemo(() => {
    if (!publicVcfFolder) {
      return { sequences: [], spanList: [], motifIdsList: [], sortedSamples: [] };
    }
    
    // Check cache first - instant return
    if (cachedData && cachedData.sequencesData.sortedSamples.length > 0) {
      return cachedData.sequencesData;
    }
    
    if (Object.keys(effectivePopulationRecords).length === 0) {
      return { sequences: [], spanList: [], motifIdsList: [], sortedSamples: [] };
    }
    
    const data = createSequencesData(record, effectivePopulationRecords);
    
    // Cache the result
    if (cacheKey) {
      const cached = processedDataCache.get(cacheKey);
      if (cached) {
        cached.sequencesData = data;
        cached.populationRecords = effectivePopulationRecords;
      } else {
        processedDataCache.set(cacheKey, {
          genotypeComparisonData: {},
          sequencesData: data,
          populationRecords: effectivePopulationRecords
        });
      }
    }
    
    return data;
  }, [record, effectivePopulationRecords, publicVcfFolder, cacheKey, cachedData]);

  // Memoize sequences data to prevent unnecessary re-renders - MUST be called before any early returns
  const memoizedSequences = useMemo(() => sequencesData.sequences, [sequencesData.sequences]);
  const memoizedSpanList = useMemo(() => sequencesData.spanList, [sequencesData.spanList]);
  const memoizedMotifIdsList = useMemo(() => sequencesData.motifIdsList, [sequencesData.motifIdsList]);
  const memoizedSortedSamples = useMemo(() => sequencesData.sortedSamples, [sequencesData.sortedSamples]);
  
  // Helper functions for sorting
  const getSampleLength = useCallback((sample: string): number => {
    const seqIdx = memoizedSequences.findIndex(s => s.name === sample);
    if (seqIdx === -1) return 0;
    return memoizedSequences[seqIdx].sequence.length;
  }, [memoizedSequences]);
  
  const getMotifCount = useCallback((sample: string): number => {
    const seqIdx = memoizedSequences.findIndex(s => s.name === sample);
    if (seqIdx === -1) return 0;
    const motifIds = memoizedMotifIdsList[seqIdx] || [];
    return motifIds.filter(id => id && id !== '.' && id !== '').length;
  }, [memoizedSequences, memoizedMotifIdsList]);
  
  const getSampleGenotype = useCallback((sample: string): string => {
    let genotype = genotypeComparisonData[sample] || '';
    if (!genotype && (sample.endsWith('_h1') || sample.endsWith('_h2'))) {
      const baseName = sample.slice(0, -3);
      genotype = genotypeComparisonData[baseName] || '';
    }
    return genotype ? String(genotype).trim() : '';
  }, [genotypeComparisonData]);
  
  // Filter samples based on search query and selected genotypes
  const filteredSamples = useMemo(() => {
    if (!searchQuery && selectedGenotypes.length === 0) {
      return memoizedSortedSamples;
    }
    
    return memoizedSortedSamples.filter(sample => {
      // Always include reference and allele samples
      if (sample === 'Ref' || sample.startsWith('Allel')) {
        const matchesSearch = !searchQuery || sample.toLowerCase().includes(searchQuery.toLowerCase());
        return matchesSearch;
      }
      
      const matchesSearch = !searchQuery || sample.toLowerCase().includes(searchQuery.toLowerCase());
      
      // Get genotype for this sample - check both direct match and potential base name match
      const genotype = getSampleGenotype(sample);
      
      // Normalize genotype for comparison (trim whitespace, convert to string)
      const normalizedSelectedGenotypes = selectedGenotypes.map(gt => String(gt).trim());
      const matchesGenotype = selectedGenotypes.length === 0 || (genotype && normalizedSelectedGenotypes.includes(genotype));
      return matchesSearch && matchesGenotype;
    });
  }, [memoizedSortedSamples, searchQuery, selectedGenotypes, getSampleGenotype]);
  
  // Sort filtered samples
  const sortedFilteredSamples = useMemo(() => {
    const samples = [...filteredSamples];
    
    // Separate reference/allele samples from population samples
    const refAlleleSamples = samples.filter(s => s === 'Ref' || s.startsWith('Allel'));
    const populationSamples = samples.filter(s => s !== 'Ref' && !s.startsWith('Allel'));
    
    // Sort population samples
    populationSamples.sort((a, b) => {
      let comparison = 0;
      switch(sortBy) {
        case 'length':
          comparison = getSampleLength(a) - getSampleLength(b);
          break;
        case 'motifCount':
          comparison = getMotifCount(a) - getMotifCount(b);
          break;
        case 'genotype':
          comparison = getSampleGenotype(a).localeCompare(getSampleGenotype(b));
          break;
        case 'name':
          comparison = a.localeCompare(b);
          break;
      }
      return sortDirection === 'asc' ? comparison : -comparison;
    });
    
    // Always keep reference/allele samples at the top
    return [...refAlleleSamples, ...populationSamples];
  }, [filteredSamples, sortBy, sortDirection, getSampleLength, getMotifCount, getSampleGenotype]);
  
  // Virtual scrolling: determine visible samples
  const shouldUseVirtualScroll = sortedFilteredSamples.length > VIRTUAL_SCROLL_THRESHOLD;
  const visibleSamples = useMemo(() => {
    if (!shouldUseVirtualScroll) {
      return sortedFilteredSamples;
    }
    
    // Always include ref/allele samples
    const refAlleleSamples = sortedFilteredSamples.filter(s => s === 'Ref' || s.startsWith('Allel'));
    const populationSamples = sortedFilteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel'));
    
    // Get visible population samples
    const visiblePopulationSamples = populationSamples.slice(visibleRange.start, visibleRange.end);
    
    return [...refAlleleSamples, ...visiblePopulationSamples];
  }, [sortedFilteredSamples, shouldUseVirtualScroll, visibleRange]);
  
  // Fetch pathogenic info for the region
  useEffect(() => {
    const fetchPathogenicInfo = async () => {
      if (!record.chr || !record.pos || !record.stop) {
        setPathogenicInfo(null);
        return;
      }
      
      try {
        const response = await axios.get(`${API_BASE}/api/pathogenic/check`, {
          params: {
            chr: record.chr,
            start: record.pos,
            end: record.stop
          }
        });
        
        console.log('Pathogenic check response:', response.data);
        if (response.data.pathogenic) {
          // Set full pathogenic info including gene, disease, inheritance
          const threshold = response.data.pathogenic_threshold;
          console.log('Pathogenic region found, threshold:', threshold);
          setPathogenicInfo({ 
            pathogenic_threshold: threshold !== undefined && threshold !== null ? threshold : undefined,
            gene: response.data.gene || undefined,
            disease: response.data.disease || undefined,
            inheritance: response.data.inheritance || undefined
          });
        } else {
          console.log('Region is not pathogenic');
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
  }, [record.chr, record.pos, record.stop]);
  
  // Update virtual scroll when samples change
  useEffect(() => {
    if (shouldUseVirtualScroll && visibleRange.end > sortedFilteredSamples.length) {
      const populationCount = sortedFilteredSamples.filter((s: string) => s !== 'Ref' && !s.startsWith('Allel')).length;
      setVisibleRange({ start: 0, end: Math.min(50, populationCount) });
    }
    setVirtualScrollEnabled(shouldUseVirtualScroll);
  }, [sortedFilteredSamples.length, shouldUseVirtualScroll]);
  
  // Filter sequences, spanList, and motifIdsList to match visible samples (for virtual scrolling)
  // Note: sequences array includes Ref, Allele1, Allele2 first, then population samples
  const filteredSequencesData = useMemo(() => {
    const samplesToUse = shouldUseVirtualScroll ? visibleSamples : sortedFilteredSamples;
    
    if (samplesToUse.length === memoizedSortedSamples.length && 
        JSON.stringify(samplesToUse) === JSON.stringify(memoizedSortedSamples)) {
      return {
        sequences: memoizedSequences,
        spanList: memoizedSpanList,
        motifIdsList: memoizedMotifIdsList
      };
    }
    
    // Create a map for quick lookup
    const sequenceMap = new Map(memoizedSequences.map((seq, idx) => [seq.name, idx]));
    
    // Get indices in the order of samplesToUse
    const orderedIndices: number[] = [];
    samplesToUse.forEach(sample => {
      const idx = sequenceMap.get(sample);
      if (idx !== undefined) {
        orderedIndices.push(idx);
      }
    });
    
    const filteredSequences = orderedIndices.map(idx => memoizedSequences[idx]);
    const filteredSpanList = orderedIndices.map(idx => memoizedSpanList[idx]);
    const filteredMotifIdsList = orderedIndices.map(idx => memoizedMotifIdsList[idx]);
    
    return {
      sequences: filteredSequences,
      spanList: filteredSpanList,
      motifIdsList: filteredMotifIdsList
    };
  }, [visibleSamples, sortedFilteredSamples, shouldUseVirtualScroll, memoizedSortedSamples, memoizedSequences, memoizedSpanList, memoizedMotifIdsList]);
  
  // Compute statistics for summary card
  const statistics = useMemo(() => {
    if (memoizedSequences.length === 0) {
      return {
        totalSamples: 0,
        uniqueGenotypes: 0,
        avgMotifCount: 0,
        minLength: 0,
        maxLength: 0,
        avgLength: 0
      };
    }
    
    const genotypeGroups: { [gt: string]: number } = {};
    Object.values(genotypeComparisonData).forEach(gt => {
      genotypeGroups[gt] = (genotypeGroups[gt] || 0) + 1;
    });
    
    const motifCounts = memoizedSequences.map((seq, idx) => {
      const motifIds = memoizedMotifIdsList[idx] || [];
      return motifIds.filter(id => id && id !== '.' && id !== '').length;
    });
    
    const lengths = memoizedSequences.map(seq => seq.sequence.length);
    
    const avgMotifCount = motifCounts.reduce((sum, count) => sum + count, 0) / motifCounts.length;
    const avgLength = lengths.reduce((sum, len) => sum + len, 0) / lengths.length;
    const minLength = Math.min(...lengths);
    const maxLength = Math.max(...lengths);
    
    return {
      totalSamples: memoizedSortedSamples.length,
      uniqueGenotypes: Object.keys(genotypeGroups).length,
      avgMotifCount: Math.round(avgMotifCount * 10) / 10,
      minLength,
      maxLength,
      avgLength: Math.round(avgLength)
    };
  }, [memoizedSequences, memoizedMotifIdsList, memoizedSortedSamples, genotypeComparisonData]);
  
  // Get unique genotypes for filter (exclude empty strings and null values)
  const availableGenotypes = useMemo(() => {
    const genotypes = Object.values(genotypeComparisonData)
      .filter(gt => gt && gt !== '' && gt !== 'null' && gt !== 'undefined')
      .map(gt => String(gt));
    return Array.from(new Set(genotypes)).sort();
  }, [genotypeComparisonData]);
  
  // Compute genotype frequency distribution
  const genotypeFrequency = useMemo(() => {
    const frequency: { [gt: string]: number } = {};
    Object.values(genotypeComparisonData).forEach(gt => {
      const normalizedGt = String(gt).trim();
      if (normalizedGt && normalizedGt !== '') {
        frequency[normalizedGt] = (frequency[normalizedGt] || 0) + 1;
      }
    });
    
    return Object.entries(frequency)
      .map(([gt, count]) => ({ genotype: gt, count, percentage: (count / Object.keys(genotypeComparisonData).length) * 100 }))
      .sort((a, b) => b.count - a.count);
  }, [genotypeComparisonData]);
  
  // Export to CSV function
  const exportToCSV = useCallback(() => {
    const headers = ['Sample', 'Genotype', 'Length (bp)', 'Motif Count', 'Sequence'];
    const rows = sortedFilteredSamples.map(sample => {
      const genotype = getSampleGenotype(sample);
      const length = getSampleLength(sample);
      const motifCount = getMotifCount(sample);
      const seqIdx = memoizedSequences.findIndex(s => s.name === sample);
      const sequence = seqIdx !== -1 ? memoizedSequences[seqIdx].sequence : '';
      
      return [sample, genotype, length.toString(), motifCount.toString(), sequence];
    });
    
    const csvContent = [
      headers.join(','),
      ...rows.map(row => row.map(cell => `"${String(cell).replace(/"/g, '""')}"`).join(','))
    ].join('\n');
    
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', `population_comparison_${region || 'data'}_${new Date().toISOString().split('T')[0]}.csv`);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
  }, [sortedFilteredSamples, getSampleGenotype, getSampleLength, getMotifCount, memoizedSequences, region]);
  
  // Filter genotype comparison data to only include filtered samples
  const filteredGenotypeComparisonData = useMemo(() => {
    if (filteredSamples.length === memoizedSortedSamples.length) {
      return genotypeComparisonData;
    }
    
    const filtered: { [sample: string]: string } = {};
    filteredSamples.forEach(sample => {
      if (genotypeComparisonData[sample] !== undefined) {
        filtered[sample] = genotypeComparisonData[sample];
      }
    });
    
    // Always include reference/allele samples
    ['Ref', 'Allel1', 'Allel2', 'Current Sample'].forEach(sample => {
      if (genotypeComparisonData[sample] !== undefined) {
        filtered[sample] = genotypeComparisonData[sample];
      }
    });
    
    return filtered;
  }, [filteredSamples, memoizedSortedSamples, genotypeComparisonData]);

  // Track last fetched region to avoid refetching
  const lastFetchedRegion = useRef<string>('');
  const lastFetchedFolder = useRef<string>('');
  
  // Check if we have any population records - MUST be called before any early returns
  const hasData = useMemo(() => {
    if (cachedData && Object.keys(cachedData.populationRecords).length > 0) {
      return true;
    }
    if (Object.keys(effectivePopulationRecords).length > 0) {
      return true;
    }
    return false;
  }, [cachedData, effectivePopulationRecords]);

  // Extract current alleles from record - MUST be called before any early returns
  const currentAlleles = useMemo(() => {
    if (!record) return [];
    const alleles: string[] = [];
    if (record.alt_allele1) alleles.push(record.alt_allele1);
    if (record.alt_allele2) alleles.push(record.alt_allele2);
    return alleles;
  }, [record]);

  // Hooks for export functionality - MUST be called before any early returns
  const populationComparisonRef = useRef<HTMLDivElement>(null);
  
  const handleExportFASTA = useCallback(() => {
    const sequences = memoizedSequences.map(seq => ({
      name: seq.name,
      sequence: seq.sequence
    }));
    exportToFASTA(sequences, generateFilename('population_sequences', 'fasta'));
  }, [memoizedSequences]);
  
  useEffect(() => {
    if (!publicVcfFolder || !region) {
      setPopulationRecords({});
      lastFetchedRegion.current = '';
      lastFetchedFolder.current = '';
      setLoading(false);
      return;
    }

    // Check cache FIRST - if we have cached data, use it immediately and skip fetch
    if (cacheKey && processedDataCache.has(cacheKey)) {
      const cached = processedDataCache.get(cacheKey)!;
      if (Object.keys(cached.populationRecords).length > 0) {
        // We have cached data - use it immediately, no loading state, no API call
        setPopulationRecords(cached.populationRecords);
        setLoading(false);
        setError(null);
        lastFetchedRegion.current = region;
        lastFetchedFolder.current = publicVcfFolder;
        return; // Skip API call entirely - instant display
      }
    }
    
    // Skip if we already have data for this exact region and folder
    if (lastFetchedRegion.current === region && lastFetchedFolder.current === publicVcfFolder) {
      return; // Already fetched, skip
    }

    // Only fetch if we don't have cached data
    const fetchPopulationData = async () => {
      setLoading(true);
      setError(null);
      try {
        const response = await axios.get(`${API_BASE}/api/population/region/${encodeURIComponent(region)}`, {
          params: { folder_path: publicVcfFolder }
        });
        const records = response.data.records || {};
        
        // Debug logging
        if (Object.keys(records).length === 0) {
          console.warn(`No population records found for region ${region} in folder ${publicVcfFolder}`);
          console.warn('Response data:', response.data);
        } else {
          console.log(`Loaded ${Object.keys(records).length} population records for region ${region}`);
        }
        
        setPopulationRecords(records);
        lastFetchedRegion.current = region;
        lastFetchedFolder.current = publicVcfFolder;
        
        // Update cache with fetched records
        if (cacheKey) {
          const cached = processedDataCache.get(cacheKey);
          if (cached) {
            cached.populationRecords = records;
          } else {
            processedDataCache.set(cacheKey, {
              genotypeComparisonData: {},
              sequencesData: { sequences: [], spanList: [], motifIdsList: [], sortedSamples: [] },
              populationRecords: records
            });
          }
        }
        
              // Pre-compute all data in parallel using Promise.all for instant subsequent renders
        // This runs asynchronously and doesn't block the UI
        if (cacheKey && record) {
          // Use requestIdleCallback if available, otherwise setTimeout for next tick
          const scheduleWork = (callback: () => void, delay: number = 0) => {
            if (typeof requestIdleCallback !== 'undefined' && delay === 0) {
              requestIdleCallback(callback, { timeout: 1000 });
            } else {
              setTimeout(callback, delay);
            }
          };
          
          // Process in parallel - compute genotype and sequences data simultaneously
          Promise.all([
            // Compute genotype comparison
            new Promise<void>((resolve) => {
              scheduleWork(() => {
                const genotypeData = computeGenotypeComparison(record, records);
                const cached = processedDataCache.get(cacheKey!);
                if (cached) {
                  cached.genotypeComparisonData = genotypeData;
                }
                resolve();
              });
            }),
            // Compute sequences data
            new Promise<void>((resolve) => {
              scheduleWork(() => {
                const seqData = createSequencesData(record, records);
                const cached = processedDataCache.get(cacheKey!);
                if (cached) {
                  cached.sequencesData = seqData;
                }
                resolve();
              });
            })
          ]).then(() => {
            // After both are done, pre-compute plot data if we have sequences
            const cached = processedDataCache.get(cacheKey!);
            if (cached && cached.sequencesData.sortedSamples.length > 0) {
              // Pre-compute plot data in background (non-blocking)
              scheduleWork(() => {
                precomputePlotData(
                  record,
                  cached.sequencesData.sequences,
                  cached.sequencesData.spanList,
                  cached.sequencesData.motifIdsList,
                  cached.sequencesData.sortedSamples
                );
              });
            }
          });
        }
      } catch (err: any) {
        console.error('Error loading population data:', err);
        console.error('Error response:', err.response?.data);
        const errorMessage = err.response?.data?.detail || err.message || 'Failed to load population data';
        setError(errorMessage);
        setPopulationRecords({});
        lastFetchedRegion.current = '';
        lastFetchedFolder.current = '';
      } finally {
        setLoading(false);
      }
    };

    fetchPopulationData();
  }, [publicVcfFolder, region, cacheKey]);

  if (!publicVcfFolder) {
    return (
      <div className="population-no-data">
        <p>No population folder specified. Please load a public VCF folder from the sidebar.</p>
      </div>
    );
  }

  if (error) {
    return <div className="population-error">Error: {error}</div>;
  }

  // Show loading only if we don't have cached data
  if (loading && !hasData) {
    return <div className="population-loading">Loading population data...</div>;
  }

  if (!hasData) {
    return (
      <div className="population-no-data">
        <p>No population data found for this region.</p>
      </div>
    );
  }

  return (
    <div className="population-comparison" ref={populationComparisonRef}>
      <div className="population-comparison-header" style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '1rem', width: '100%', gap: '1rem' }}>
        <h2 style={{ margin: 0, flex: '1 1 auto', minWidth: 0 }}>🌐 Population Comparison</h2>
        <div style={{ flexShrink: 0 }}>
          <ExportMenu
            elementRef={populationComparisonRef}
            filename="population_comparison"
            onExportCSV={exportToCSV}
            onExportFASTA={handleExportFASTA}
            showImageExport={true}
            showDataExport={true}
          />
        </div>
      </div>
      
      {record && (
        <PopulationFrequencyPanel
          chr={record.chr}
          pos={record.pos}
          stop={record.stop}
          region={region}
          currentAlleles={currentAlleles}
        />
      )}
      
      {/* Enhanced Statistics Summary Card */}
      <div className="population-stats-summary-enhanced">
        <div className="stats-header">
          <h3 className="stats-title">📈 Population Statistics</h3>
        </div>
        <div className="stats-grid">
          <div className="stat-card-enhanced">
            <div className="stat-icon">👥</div>
            <div className="stat-content">
              <span className="stat-label-enhanced">Total Samples</span>
              <span className="stat-value-enhanced">{statistics.totalSamples}</span>
            </div>
          </div>
          <div className="stat-card-enhanced">
            <div className="stat-icon">🧬</div>
            <div className="stat-content">
              <span className="stat-label-enhanced">Unique Genotypes</span>
              <span className="stat-value-enhanced">{statistics.uniqueGenotypes}</span>
            </div>
          </div>
          <div className="stat-card-enhanced">
            <div className="stat-icon">🔢</div>
            <div className="stat-content">
              <span className="stat-label-enhanced">Avg Motif Count</span>
              <span className="stat-value-enhanced">{statistics.avgMotifCount}</span>
            </div>
          </div>
          <div className="stat-card-enhanced">
            <div className="stat-icon">📏</div>
            <div className="stat-content">
              <span className="stat-label-enhanced">Length Range</span>
              <span className="stat-value-enhanced">{statistics.minLength}-{statistics.maxLength}bp</span>
            </div>
          </div>
          <div className="stat-card-enhanced">
            <div className="stat-icon">📊</div>
            <div className="stat-content">
              <span className="stat-label-enhanced">Avg Length</span>
              <span className="stat-value-enhanced">{statistics.avgLength}bp</span>
            </div>
          </div>
        </div>
      </div>
      
      {/* Genotype Frequency Distribution */}
      {/* Pathogenic Information Panel */}
      {pathogenicInfo && (pathogenicInfo.gene || pathogenicInfo.disease || pathogenicInfo.inheritance) && (
        <div className="population-section pathogenic-info-section">
          <div className="pathogenic-info-header">
            <h3>⚠️ Pathogenic Region Information</h3>
          </div>
          <div className="pathogenic-info-content">
            {pathogenicInfo.gene && (
              <div className="pathogenic-info-item">
                <span className="pathogenic-info-label">Gene:</span>
                <span className="pathogenic-info-value">{pathogenicInfo.gene}</span>
              </div>
            )}
            {pathogenicInfo.disease && (
              <div className="pathogenic-info-item">
                <span className="pathogenic-info-label">Disease:</span>
                <span className="pathogenic-info-value">{pathogenicInfo.disease}</span>
              </div>
            )}
            {pathogenicInfo.inheritance && (
              <div className="pathogenic-info-item">
                <span className="pathogenic-info-label">Inheritance:</span>
                <span className="pathogenic-info-value">{pathogenicInfo.inheritance}</span>
              </div>
            )}
            {pathogenicInfo.pathogenic_threshold && (
              <div className="pathogenic-info-item">
                <span className="pathogenic-info-label">Pathogenic Threshold:</span>
                <span className="pathogenic-info-value threshold">{pathogenicInfo.pathogenic_threshold} repeats</span>
              </div>
            )}
          </div>
        </div>
      )}

      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, frequency: !expandedSections.frequency})}
        >
          <h3>Genotype Frequency Distribution</h3>
          <span className="collapse-icon">{expandedSections.frequency ? '▼' : '▶'}</span>
        </div>
        {expandedSections.frequency && (
          <div className="frequency-chart">
            <div className="frequency-bars">
              {genotypeFrequency.map(({ genotype, count, percentage }) => (
                <div key={genotype} className="frequency-bar-item">
                  <div className="frequency-bar-label">
                    <span className="genotype-label">{genotype}</span>
                    <span className="frequency-count">{count} ({percentage.toFixed(1)}%)</span>
                  </div>
                  <div className="frequency-bar-container">
                    <div 
                      className="frequency-bar"
                      style={{ width: `${percentage}%` }}
                      title={`${genotype}: ${count} samples (${percentage.toFixed(1)}%)`}
                    />
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Genotype Comparison Matrix */}
      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, genotype: !expandedSections.genotype})}
        >
          <h3>Genotype Comparison Matrix</h3>
          <span className="collapse-icon">{expandedSections.genotype ? '▼' : '▶'}</span>
        </div>
        {expandedSections.genotype && (
          <GenotypeComparisonMatrix genotypes={filteredGenotypeComparisonData} />
        )}
      </div>

      {/* Search and Filter Bar - Above Stack Plot and Heatmap */}
      <div className="population-search-filter">
        <div className="search-wrapper">
          <input
            type="text"
            className="population-search-input"
            placeholder="🔍 Search samples..."
            value={searchQuery}
            onChange={(e) => setSearchQuery(e.target.value)}
          />
        </div>
        <div className="filter-wrapper">
          <div className="filter-label">Filter by Genotype:</div>
          <div className="genotype-filter-buttons">
            {availableGenotypes.map(gt => (
              <button
                key={gt}
                className={`genotype-filter-btn ${selectedGenotypes.includes(gt) ? 'active' : ''}`}
                onClick={() => {
                  if (selectedGenotypes.includes(gt)) {
                    setSelectedGenotypes(selectedGenotypes.filter(g => g !== gt));
                  } else {
                    setSelectedGenotypes([...selectedGenotypes, gt]);
                  }
                }}
              >
                {gt}
              </button>
            ))}
            {selectedGenotypes.length > 0 && (
              <button
                className="genotype-filter-btn clear"
                onClick={() => setSelectedGenotypes([])}
              >
                Clear
              </button>
            )}
          </div>
        </div>
        {sortedFilteredSamples.length !== memoizedSortedSamples.length && (
          <div className="filter-indicator">
            Showing {sortedFilteredSamples.length} of {memoizedSortedSamples.length} samples
          </div>
        )}
        
        {/* Sort Controls */}
        <div className="sort-controls">
          <div className="sort-label">Sort by:</div>
          <div className="sort-buttons">
            <button
              className={`sort-btn ${sortBy === 'length' ? 'active' : ''}`}
              onClick={() => {
                if (sortBy === 'length') {
                  setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                } else {
                  setSortBy('length');
                  setSortDirection('desc');
                }
              }}
            >
              Length {sortBy === 'length' && (sortDirection === 'asc' ? '↑' : '↓')}
            </button>
            <button
              className={`sort-btn ${sortBy === 'motifCount' ? 'active' : ''}`}
              onClick={() => {
                if (sortBy === 'motifCount') {
                  setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                } else {
                  setSortBy('motifCount');
                  setSortDirection('desc');
                }
              }}
            >
              Motif Count {sortBy === 'motifCount' && (sortDirection === 'asc' ? '↑' : '↓')}
            </button>
            <button
              className={`sort-btn ${sortBy === 'name' ? 'active' : ''}`}
              onClick={() => {
                if (sortBy === 'name') {
                  setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                } else {
                  setSortBy('name');
                  setSortDirection('asc');
                }
              }}
            >
              Name {sortBy === 'name' && (sortDirection === 'asc' ? '↑' : '↓')}
            </button>
            <button
              className={`sort-btn ${sortBy === 'genotype' ? 'active' : ''}`}
              onClick={() => {
                if (sortBy === 'genotype') {
                  setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                } else {
                  setSortBy('genotype');
                  setSortDirection('asc');
                }
              }}
            >
              Genotype {sortBy === 'genotype' && (sortDirection === 'asc' ? '↑' : '↓')}
            </button>
          </div>
        </div>
      </div>

      {/* Stack Plot and Heatmap */}
      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, stackHeatmap: !expandedSections.stackHeatmap})}
        >
          <h3>Motif Stack Plot & Heatmap</h3>
          <span className="collapse-icon">{expandedSections.stackHeatmap ? '▼' : '▶'}</span>
        </div>
        {expandedSections.stackHeatmap && (
          <>
            {virtualScrollEnabled && (
              <div className="virtual-scroll-info">
                Showing {visibleSamples.length} of {sortedFilteredSamples.length} samples
                <div className="virtual-scroll-controls">
                  <button
                    className="virtual-scroll-btn"
                    onClick={() => {
                      const populationCount = sortedFilteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length;
                      setVisibleRange(prev => ({
                        start: Math.max(0, prev.start - 20),
                        end: Math.min(populationCount, prev.end - 20)
                      }));
                    }}
                    disabled={visibleRange.start === 0}
                  >
                    ↑ Load More Above
                  </button>
                  <button
                    className="virtual-scroll-btn"
                    onClick={() => {
                      const populationCount = sortedFilteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length;
                      setVisibleRange(prev => ({
                        start: prev.start + 20,
                        end: Math.min(populationCount, prev.end + 20)
                      }));
                    }}
                    disabled={visibleRange.end >= sortedFilteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length}
                  >
                    ↓ Load More Below
                  </button>
                </div>
              </div>
            )}
            <CombinedStackHeatmap 
              record={record}
              sequences={filteredSequencesData.sequences}
              spanList={filteredSequencesData.spanList}
              motifIdsList={filteredSequencesData.motifIdsList}
              sortedSamples={visibleSamples}
              pathogenicThreshold={pathogenicInfo?.pathogenic_threshold}
            />
          </>
        )}
      </div>

      {/* Bar Plot */}
      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, barPlot: !expandedSections.barPlot})}
        >
          <h3>Motif Count per Sample</h3>
          <span className="collapse-icon">{expandedSections.barPlot ? '▼' : '▶'}</span>
        </div>
        {expandedSections.barPlot && (
          <BarPlot 
            sequences={filteredSequencesData.sequences} 
            motifIdsList={filteredSequencesData.motifIdsList}
            pathogenicThreshold={pathogenicInfo?.pathogenic_threshold}
          />
        )}
      </div>

      {/* Cluster Plot */}
      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, cluster: !expandedSections.cluster})}
        >
          <h3>Sample Clustering</h3>
          <span className="collapse-icon">{expandedSections.cluster ? '▼' : '▶'}</span>
        </div>
        {expandedSections.cluster && (
          <ClusterPlot sequences={filteredSequencesData.sequences} motifIdsList={filteredSequencesData.motifIdsList} />
        )}
      </div>
    </div>
  );
};

// Helper function to compute genotype comparison
function computeGenotypeComparison(
  record: Record,
  populationRecords: { [key: string]: PopulationRecord }
): { [sample: string]: string } {
  const genotypes: { [sample: string]: string } = {};
  
  // Add current sample
  genotypes['Current Sample'] = record.gt;
  
  // Group h1 and h2 records by base sample name
  const sampleGroups: { [base: string]: { h1?: PopulationRecord; h2?: PopulationRecord } } = {};
  
  for (const [sampleName, popRecord] of Object.entries(populationRecords)) {
    if (sampleName.endsWith('_h1')) {
      const baseName = sampleName.slice(0, -3);
      if (!sampleGroups[baseName]) sampleGroups[baseName] = {};
      sampleGroups[baseName].h1 = popRecord;
    } else if (sampleName.endsWith('_h2')) {
      const baseName = sampleName.slice(0, -3);
      if (!sampleGroups[baseName]) sampleGroups[baseName] = {};
      sampleGroups[baseName].h2 = popRecord;
    } else {
      // Single haplotype file
      genotypes[sampleName] = popRecord.gt || '0';
    }
  }
  
  // Compute diploid genotypes from h1/h2 pairs
  for (const [baseName, group] of Object.entries(sampleGroups)) {
    if (group.h1 && group.h2) {
      // Both haplotypes available - compute diploid genotype
      const gt_h1 = group.h1.gt || '0';
      const gt_h2 = group.h2.gt || '0';
      const ids_h1 = group.h1.motif_ids_h || [];
      const ids_h2 = group.h2.motif_ids_h || [];
      
      genotypes[baseName] = computeDiploidGenotype(gt_h1, gt_h2, ids_h1, ids_h2);
    } else if (group.h1) {
      genotypes[baseName] = group.h1.gt || '0';
    } else if (group.h2) {
      genotypes[baseName] = group.h2.gt || '0';
    }
  }
  
  return genotypes;
}

// Compute diploid genotype from two haplotypes
function computeDiploidGenotype(gt_h1: string, gt_h2: string, ids_h1: string[], ids_h2: string[]): string {
  const g1 = parseInt(gt_h1) || 0;
  const g2 = parseInt(gt_h2) || 0;
  
  if (g1 === 0 && g2 === 0) return '0/0';
  if ((g1 === 0 && g2 === 1) || (g1 === 1 && g2 === 0)) return '0/1';
  if (g1 === 1 && g2 === 1) {
    // Compare motif IDs to determine if homozygous (1/1) or heterozygous (1/2)
    const ids1Str = ids_h1.join('_');
    const ids2Str = ids_h2.join('_');
    return ids1Str === ids2Str ? '1/1' : '1/2';
  }
  return `${g1}/${g2}`;
}

// Create sequences data for plots - optimized with chunked processing
function createSequencesData(
  record: Record,
  populationRecords: { [key: string]: PopulationRecord }
) {
  const sequences: Array<{ name: string; sequence: string }> = [];
  const spanList: string[] = [];
  const motifIdsList: string[][] = [];
  
  // Add reference
  sequences.push({ name: 'Ref', sequence: record.ref_allele });
  spanList.push(record.spans[0] || '');
  motifIdsList.push(record.motif_ids_ref || []);
  
  // Add allele 1
  sequences.push({ name: 'Allel1', sequence: record.alt_allele1 });
  spanList.push(record.spans[1] || '');
  motifIdsList.push(record.motif_ids_h1 || []);
  
  // Add allele 2 if present
  if (record.alt_allele2 && record.alt_allele2 !== '') {
    sequences.push({ name: 'Allel2', sequence: record.alt_allele2 });
    spanList.push(record.spans[2] || '');
    motifIdsList.push(record.motif_ids_h2 || []);
  }
  
  // Add population samples - optimized with direct iteration
  const popEntries = Object.entries(populationRecords);
  for (let i = 0; i < popEntries.length; i++) {
    const [sampleName, popRecord] = popEntries[i];
    if (popRecord.alt_allele && popRecord.alt_allele !== '') {
      sequences.push({ name: sampleName, sequence: popRecord.alt_allele });
      spanList.push(popRecord.spans || '');
      motifIdsList.push(popRecord.motif_ids_h || []);
    }
  }
  
  // Calculate sorted samples by total length (for consistent ordering across plots)
  // Optimize: pre-compute lengths with single-pass calculation
  // Use chunked processing for large datasets to avoid blocking
  const sampleLengths: Array<{ sample: string; totalLength: number }> = [];
  const CHUNK_SIZE = 50; // Process 50 samples at a time
  
  // Process in chunks if we have many samples
  if (sequences.length > CHUNK_SIZE) {
    // For large datasets, we'll process synchronously but in optimized chunks
    // This is still faster than the original approach
    for (let chunkStart = 0; chunkStart < sequences.length; chunkStart += CHUNK_SIZE) {
      const chunkEnd = Math.min(chunkStart + CHUNK_SIZE, sequences.length);
      for (let idx = chunkStart; idx < chunkEnd; idx++) {
        const span = spanList[idx] || '';
        const ranges = parseMotifRange(span); // Now cached
        // Optimize: single-pass sum calculation
        let totalLength = 0;
        for (let r = 0; r < ranges.length; r++) {
          const [start, end] = ranges[r];
          totalLength += (end + 1 - start);
        }
        sampleLengths.push({ sample: sequences[idx].name, totalLength });
      }
    }
  } else {
    // For smaller datasets, process all at once
    for (let idx = 0; idx < sequences.length; idx++) {
      const span = spanList[idx] || '';
      const ranges = parseMotifRange(span); // Now cached
      // Optimize: single-pass sum calculation
      let totalLength = 0;
      for (let r = 0; r < ranges.length; r++) {
        const [start, end] = ranges[r];
        totalLength += (end + 1 - start);
      }
      sampleLengths.push({ sample: sequences[idx].name, totalLength });
    }
  }
  
  // Sort in-place for better performance
  sampleLengths.sort((a, b) => b.totalLength - a.totalLength);
  const sortedSamples = sampleLengths.map(s => s.sample);
  
  return { sequences, spanList, motifIdsList, sortedSamples };
}

// Pre-compute plot data asynchronously to avoid blocking render
function precomputePlotData(
  record: Record,
  sequences: Array<{ name: string; sequence: string }>,
  spanList: string[],
  motifIdsList: string[][],
  sortedSamples: string[]
) {
  // Create cache key
  const seqKey = sequences.map(s => s.name).join(',');
  const spanKey = spanList.join('|');
  const motifKey = motifIdsList.map(m => m.join('_')).join('|');
  const plotCacheKey = `${record.id}:${seqKey}:${spanKey}:${motifKey}`;
  
  // Check if already cached
  if (combinedPlotCache.has(plotCacheKey)) {
    return; // Already computed
  }
  
  // Process motif data for stack plot - optimized with chunked processing
  const stackData: Array<{
    sample: string;
    start: number;
    end: number;
    motif: string;
    motifIndex: number;
    length: number;
  }> = [];
  
  const CHUNK_SIZE = 20; // Process 20 sequences at a time
  
  // Process sequences in chunks
  for (let chunkStart = 0; chunkStart < sequences.length; chunkStart += CHUNK_SIZE) {
    const chunkEnd = Math.min(chunkStart + CHUNK_SIZE, sequences.length);
    
    for (let seqIdx = chunkStart; seqIdx < chunkEnd; seqIdx++) {
      const seq = sequences[seqIdx];
      const span = spanList[seqIdx] || '';
      const motifIds = motifIdsList[seqIdx] || [];
      const ranges = parseMotifRange(span);
      
      let previousEnd = 0;
      for (let rangeIdx = 0; rangeIdx < ranges.length; rangeIdx++) {
        const [start, end] = ranges[rangeIdx];
        
        if (start > previousEnd) {
          stackData.push({
            sample: seq.name,
            start: previousEnd,
            end: start,
            motif: 'Interruption',
            motifIndex: -1,
            length: start - previousEnd
          });
        }
        
        const motifId = parseInt(motifIds[rangeIdx] || '0');
        const motifName = record.motifs[motifId] || 'Unknown';
        stackData.push({
          sample: seq.name,
          start: start,
          end: end + 1,
          motif: motifName,
          motifIndex: motifId,
          length: end + 1 - start
        });
        
        previousEnd = end + 1;
      }
      
      if (previousEnd < seq.sequence.length) {
        stackData.push({
          sample: seq.name,
          start: previousEnd,
          end: seq.sequence.length,
          motif: 'Interruption',
          motifIndex: -1,
          length: seq.sequence.length - previousEnd
        });
      }
    }
  }
  
  // Process heatmap data - optimized with pre-allocated arrays
  const heatData: Array<{ sample: string; motif: string; count: number }> = [];
  const heatDataMap = new Map<string, number>();
  
  for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
    const seq = sequences[seqIdx];
    const motifIds = motifIdsList[seqIdx] || [];
    const motifCounts: { [motif: string]: number } = {};
    
    for (let i = 0; i < motifIds.length; i++) {
      const motifId = motifIds[i];
      if (motifId && motifId !== '.' && motifId !== '') {
        const motifIdx = parseInt(motifId);
        if (!isNaN(motifIdx) && motifIdx >= 0 && motifIdx < record.motifs.length) {
          const motifName = record.motifs[motifIdx];
          motifCounts[motifName] = (motifCounts[motifName] || 0) + 1;
        }
      }
    }
    
    for (const [motif, count] of Object.entries(motifCounts)) {
      heatData.push({ sample: seq.name, motif, count });
    }
  }

  // Get unique motifs and sort - optimized
  const motifTotalsMap = new Map<string, number>();
  for (const d of heatData) {
    motifTotalsMap.set(d.motif, (motifTotalsMap.get(d.motif) || 0) + d.count);
  }
  
  const motifTotals = Array.from(motifTotalsMap.entries())
    .map(([motif, total]) => ({ motif, total }))
    .sort((a, b) => b.total - a.total);
  const sortedMots = motifTotals.map(m => m.motif);
  
  // Color mapping - optimized: use Set directly and single pass
  const uniqueStackMotifsSet = new Set<string>();
  for (const d of stackData) {
    if (d.motif !== 'Interruption') {
      uniqueStackMotifsSet.add(d.motif);
    }
  }
  const uniqueStackMotifs = Array.from(uniqueStackMotifsSet);
  const colorMap: { [key: string]: string } = {};
  for (let idx = 0; idx < uniqueStackMotifs.length; idx++) {
    const motif = uniqueStackMotifs[idx];
    const motifIdx = record.motifs.indexOf(motif);
    if (motifIdx >= 0) {
      colorMap[motif] = COLOR_PALETTE[motifIdx % COLOR_PALETTE.length];
    } else {
      colorMap[motif] = COLOR_PALETTE[idx % COLOR_PALETTE.length];
    }
  }
  colorMap['Interruption'] = '#ef4444';
  
  // Optimize: calculate max in single pass
  let maxLen = 1;
  let maxCnt = 1;
  for (const d of stackData) {
    if (d.end > maxLen) maxLen = d.end;
  }
  for (const d of heatData) {
    if (d.count > maxCnt) maxCnt = d.count;
  }
  
  // Pre-compute sample data for faster rendering
  const sampleDataMap: { [sample: string]: Array<{
    start: number;
    end: number;
    motif: string;
    motifIndex: number;
    length: number;
  }> } = {};
  
  // Initialize all samples with empty arrays
  for (const sample of sortedSamples) {
    sampleDataMap[sample] = [];
  }
  
  // Build map in one pass
  for (const item of stackData) {
    if (sampleDataMap[item.sample]) {
      sampleDataMap[item.sample].push(item);
    }
  }
  
  // Sort each sample's data once
  for (const sample of sortedSamples) {
    sampleDataMap[sample].sort((a, b) => a.start - b.start);
  }
  
  const result = {
    motifData: stackData,
    heatmapData: heatData,
    motifColorMap: colorMap,
    maxLength: maxLen,
    maxCount: maxCnt,
    sortedMotifs: sortedMots,
    sampleDataMap: sampleDataMap
  };
  
  // Cache the result
  combinedPlotCache.set(plotCacheKey, result);
  
  // Limit cache size to prevent memory issues (keep last 50 entries)
  if (combinedPlotCache.size > 50) {
    const firstKey = combinedPlotCache.keys().next().value;
    if (firstKey) {
      combinedPlotCache.delete(firstKey);
    }
  }
}

// Genotype Comparison Matrix Component
interface GenotypeComparisonMatrixProps {
  genotypes: { [sample: string]: string };
}

const GenotypeComparisonMatrix: React.FC<GenotypeComparisonMatrixProps> = memo(({ genotypes }) => {
  // Group genotypes
  const genotypeGroups: { [gt: string]: string[] } = {};
  for (const [sample, gt] of Object.entries(genotypes)) {
    if (!genotypeGroups[gt]) {
      genotypeGroups[gt] = [];
    }
    genotypeGroups[gt].push(sample);
  }
  
  // Sort groups by genotype
  const sortedGroups = Object.entries(genotypeGroups).sort(([a], [b]) => a.localeCompare(b));
  
  // Get genotype info
  const getGenotypeInfo = (gt: string) => {
    if (gt === '0/0') return { icon: '⚪', description: 'Homozygous Reference', color: '#10b981', bgColor: '#d1fae5' };
    if (gt === '0/1') return { icon: '🟡', description: 'Heterozygous', color: '#f59e0b', bgColor: '#fef3c7' };
    if (gt === '1/1') return { icon: '🟠', description: 'Homozygous Alternate', color: '#f97316', bgColor: '#ffedd5' };
    if (gt === '1/2') return { icon: '🔴', description: 'Heterozygous Alternate', color: '#ef4444', bgColor: '#fee2e2' };
    return { icon: '❓', description: 'Unknown', color: '#6b7280', bgColor: '#f3f4f6' };
  };
  
  // Calculate summary stats
  const totalSamples = Object.keys(genotypes).length;
  const uniqueGenotypes = Object.keys(genotypeGroups).length;
  const homozygous = Object.keys(genotypeGroups).filter(gt => gt === '0/0' || gt === '1/1').length;
  const heterozygous = Object.keys(genotypeGroups).filter(gt => gt === '0/1' || gt === '1/2').length;
  const missing = 0; // Could calculate if needed
  
  return (
    <div className="genotype-matrix">
      {/* Summary Statistics */}
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
      
      {/* Grouped View */}
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

// Helper function to parse motif ranges from span string
function parseMotifRange(spanStr: string): Array<[number, number]> {
  if (!spanStr) return [];
  
  // Check cache first
  if (parsedRangeCache.has(spanStr)) {
    return parsedRangeCache.get(spanStr)!;
  }
  
  const pattern = /\((\d+)-(\d+)\)/g;
  const ranges: Array<[number, number]> = [];
  let match: RegExpExecArray | null;
  
  // Use exec in a loop instead of matchAll for better compatibility
  while ((match = pattern.exec(spanStr)) !== null) {
    ranges.push([parseInt(match[1]) - 1, parseInt(match[2]) - 1]);
  }
  
  // Cache the result (limit cache size to prevent memory issues)
  if (parsedRangeCache.size > 1000) {
    const firstKey = parsedRangeCache.keys().next().value;
    if (firstKey) parsedRangeCache.delete(firstKey);
  }
  parsedRangeCache.set(spanStr, ranges);
  
  return ranges;
}

// Stack Plot Component
interface StackPlotProps {
  record: Record;
  sequences: Array<{ name: string; sequence: string }>;
  spanList: string[];
  motifIdsList: string[][];
  sortedSamples?: string[];
  pathogenicThreshold?: number; // Pathogenic threshold in copy number
}

const StackPlot: React.FC<StackPlotProps> = memo(({ record, sequences, spanList, motifIdsList, sortedSamples: providedSortedSamples, pathogenicThreshold }) => {
  // Create motif data for each sequence
  const motifData = useMemo(() => {
    const data: Array<{
      sample: string;
      start: number;
      end: number;
      motif: string;
      motifIndex: number;
      length: number;
    }> = [];
    
    sequences.forEach((seq, seqIdx) => {
      const span = spanList[seqIdx] || '';
      const motifIds = motifIdsList[seqIdx] || [];
      const ranges = parseMotifRange(span);
      
      let previousEnd = 0;
      ranges.forEach(([start, end], rangeIdx) => {
        // Add interruption if there's a gap
        if (start > previousEnd) {
          data.push({
            sample: seq.name,
            start: previousEnd,
            end: start,
            motif: 'Interruption',
            motifIndex: -1,
            length: start - previousEnd
          });
        }
        
        // Add motif segment
        const motifId = parseInt(motifIds[rangeIdx] || '0');
        const motifName = record.motifs[motifId] || 'Unknown';
        data.push({
          sample: seq.name,
          start: start,
          end: end + 1,
          motif: motifName,
          motifIndex: motifId,
          length: end + 1 - start
        });
        
        previousEnd = end + 1;
      });
      
      // Add trailing interruption if any
      if (previousEnd < seq.sequence.length) {
        data.push({
          sample: seq.name,
          start: previousEnd,
          end: seq.sequence.length,
          motif: 'Interruption',
          motifIndex: -1,
          length: seq.sequence.length - previousEnd
        });
      }
    });
    
    return data;
  }, [sequences, spanList, motifIdsList, record.motifs]);
  
  // Use provided sorted samples if available, otherwise calculate
  const sortedSamples = useMemo(() => {
    if (providedSortedSamples && providedSortedSamples.length > 0) {
      return providedSortedSamples;
    }
    // Fallback: calculate sorting
    const uniqueSamples = Array.from(new Set(motifData.map(d => d.sample)));
    const sampleLengths = uniqueSamples.map(sample => {
      const sampleData = motifData.filter(d => d.sample === sample && d.motif !== 'Interruption');
      return {
        sample,
        totalLength: sampleData.reduce((sum, d) => sum + d.length, 0)
      };
    }).sort((a, b) => b.totalLength - a.totalLength);
    return sampleLengths.map(s => s.sample);
  }, [providedSortedSamples, motifData]);
  
  // Get unique motifs for color mapping
  const uniqueMotifs = Array.from(new Set(motifData.filter(d => d.motif !== 'Interruption').map(d => d.motif)));
  const motifColorMap: { [key: string]: string } = {};
  uniqueMotifs.forEach((motif, idx) => {
    const motifIdx = record.motifs.indexOf(motif);
    if (motifIdx >= 0) {
      motifColorMap[motif] = COLOR_PALETTE[motifIdx % COLOR_PALETTE.length];
    } else {
      motifColorMap[motif] = COLOR_PALETTE[idx % COLOR_PALETTE.length];
    }
  });
  motifColorMap['Interruption'] = '#e5e7eb';
  
  // Get max length for scaling
  const maxLength = Math.max(...motifData.map(d => d.end), 1);
  
  // Calculate pathogenic threshold length (if threshold is in copy number, convert to length)
  // Assuming average motif length for conversion
  const avgMotifLength = record.motifs.length > 0 ? record.motifs[0].length : 1;
  const pathogenicThresholdLength = pathogenicThreshold && pathogenicThreshold > 0 
    ? pathogenicThreshold * avgMotifLength 
    : null;
  
  // Find samples that exceed the threshold
  const samplesExceedingThreshold = useMemo(() => {
    if (!pathogenicThresholdLength) return new Set<string>();
    const exceeding = new Set<string>();
    sortedSamples.forEach(sample => {
      const sampleData = motifData.filter(d => d.sample === sample && d.motif !== 'Interruption');
      const totalLength = sampleData.reduce((sum, d) => sum + d.length, 0);
      if (totalLength > pathogenicThresholdLength) {
        exceeding.add(sample);
      }
    });
    return exceeding;
  }, [sortedSamples, motifData, pathogenicThresholdLength]);
  
  return (
    <div className="stack-plot-container">
      <div className="stack-plot-chart" style={{ position: 'relative' }}>
        {sortedSamples.map(sample => {
          const sampleData = motifData.filter(d => d.sample === sample).sort((a, b) => a.start - b.start);
          const totalLength = Math.max(...sampleData.map(d => d.end), 0);
          const exceedsThreshold = samplesExceedingThreshold.has(sample);
          
          return (
            <div key={sample} className="stack-plot-row" style={{ position: 'relative' }}>
              <div className="stack-plot-sample-label">{sample}</div>
              <div className="stack-plot-bars-container" style={{ position: 'relative' }}>
                {/* Pathogenic threshold line */}
                {pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0 && (() => {
                  const thresholdNum = typeof pathogenicThreshold === 'number' ? pathogenicThreshold : parseFloat(String(pathogenicThreshold));
                  if (isNaN(thresholdNum) || thresholdNum <= 0) return null;
                  
                  const avgMotifLength = record.motifs.length > 0 ? record.motifs[0].length : 1;
                  const pathogenicThresholdLength = thresholdNum * avgMotifLength;
                  const isOutsideRange = pathogenicThresholdLength > maxLength;
                  const thresholdPercent = isOutsideRange 
                    ? 100 // Show at right edge if outside range
                    : (pathogenicThresholdLength / maxLength) * 100;
                  
                  return (
                    <div
                      className="pathogenic-threshold-line"
                      style={{
                        position: 'absolute',
                        left: `${thresholdPercent}%`,
                        top: 0,
                        bottom: 0,
                        width: '2px',
                        backgroundColor: '#EF4444',
                        borderLeft: isOutsideRange ? '2px solid #EF4444' : '2px dashed #EF4444',
                        zIndex: 10,
                        pointerEvents: 'none',
                        opacity: isOutsideRange ? 0.5 : 1
                      }}
                      title={`Pathogenic Threshold: ${thresholdNum} copies (${pathogenicThresholdLength}bp)${isOutsideRange ? ' (beyond visible range)' : ''}`}
                    />
                  );
                })()}
                <div className="stack-plot-bars" style={{ width: `${(totalLength / maxLength) * 100}%` }}>
                  {sampleData.map((segment, idx) => {
                    const widthPercent = (segment.length / totalLength) * 100;
                    const color = motifColorMap[segment.motif] || '#9ca3af';
                    return (
                      <div
                        key={idx}
                        className="stack-plot-segment"
                        style={{
                          width: `${widthPercent}%`,
                          backgroundColor: color,
                          borderRight: idx < sampleData.length - 1 ? '1px solid rgba(0,0,0,0.1)' : 'none'
                        }}
                        title={`${segment.motif}: ${segment.start}-${segment.end} (${segment.length}bp)`}
                      />
                    );
                  })}
                </div>
                {/* Pathogenic warning indicator */}
                {exceedsThreshold && (
                  <div
                    style={{
                      position: 'absolute',
                      left: `${(totalLength / maxLength) * 100}%`,
                      top: '50%',
                      transform: 'translateY(-50%)',
                      marginLeft: '8px',
                      color: '#DC2626',
                      fontWeight: 'bold',
                      fontSize: '14px',
                      whiteSpace: 'nowrap',
                      zIndex: 11
                    }}
                    title="Exceeds pathogenic threshold"
                  >
                    🚨 PATHOGENIC
                  </div>
                )}
              </div>
              <div className="stack-plot-length-label">{totalLength}bp</div>
            </div>
          );
        })}
      </div>
      
      {/* Legend */}
      <div className="stack-plot-legend">
        {uniqueMotifs.map(motif => (
          <div key={motif} className="legend-item">
            <span 
              className="legend-color" 
              style={{ backgroundColor: motifColorMap[motif] }}
            />
            <span className="legend-text">{motif}</span>
          </div>
        ))}
        <div className="legend-item">
          <span className="legend-color" style={{ backgroundColor: '#e5e7eb' }} />
          <span className="legend-text">Interruption</span>
        </div>
      </div>
    </div>
  );
});
StackPlot.displayName = 'StackPlot';

// Combined Stack Plot and Heatmap Component (for perfect alignment)
interface CombinedStackHeatmapProps {
  record: Record;
  sequences: Array<{ name: string; sequence: string }>;
  spanList: string[];
  motifIdsList: string[][];
  sortedSamples: string[];
  pathogenicThreshold?: number;
}

const CombinedStackHeatmap: React.FC<CombinedStackHeatmapProps> = memo(({ 
  record, 
  sequences, 
  spanList, 
  motifIdsList, 
  sortedSamples,
  pathogenicThreshold
}) => {
  // Create cache key based on inputs - MUST be called before any early returns
  const plotCacheKey = useMemo(() => {
    const seqKey = sequences.map(s => s.name).join(',');
    const spanKey = spanList.join('|');
    const motifKey = motifIdsList.map(m => m.join('_')).join('|');
    return `${record.id}:${seqKey}:${spanKey}:${motifKey}`;
  }, [record.id, sequences, spanList, motifIdsList]);
  
  // Pre-compute all data once for performance, with caching - MUST be called before any early returns
  const plotData = useMemo(() => {
    // Check cache first - instant return if cached
    if (plotCacheKey) {
      const cached = combinedPlotCache.get(plotCacheKey);
      if (cached) {
        return cached; // Instant return from cache
      }
    }
    // Process motif data for stack plot
    const stackData: Array<{
      sample: string;
      start: number;
      end: number;
      motif: string;
      motifIndex: number;
      length: number;
    }> = [];
    
    // Optimize: use for loops instead of forEach for better performance
    for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
      const seq = sequences[seqIdx];
      const span = spanList[seqIdx] || '';
      const motifIds = motifIdsList[seqIdx] || [];
      const ranges = parseMotifRange(span);
      
      // If no ranges and sequence is empty, add a single empty segment to ensure sample appears
      if (ranges.length === 0 && (!seq.sequence || seq.sequence.length === 0)) {
        stackData.push({
          sample: seq.name,
          start: 0,
          end: 0,
          motif: 'Interruption',
          motifIndex: -1,
          length: 0
        });
        continue;
      }
      
      let previousEnd = 0;
      for (let rangeIdx = 0; rangeIdx < ranges.length; rangeIdx++) {
        const [start, end] = ranges[rangeIdx];
        
        if (start > previousEnd) {
          stackData.push({
            sample: seq.name,
            start: previousEnd,
            end: start,
            motif: 'Interruption',
            motifIndex: -1,
            length: start - previousEnd
          });
        }
        
        const motifId = parseInt(motifIds[rangeIdx] || '0');
        const motifName = record.motifs[motifId] || 'Unknown';
        stackData.push({
          sample: seq.name,
          start: start,
          end: end + 1,
          motif: motifName,
          motifIndex: motifId,
          length: end + 1 - start
        });
        
        previousEnd = end + 1;
      }
      
      // Always add trailing segment if sequence has length
      if (seq.sequence && seq.sequence.length > 0) {
        if (previousEnd < seq.sequence.length) {
          stackData.push({
            sample: seq.name,
            start: previousEnd,
            end: seq.sequence.length,
            motif: 'Interruption',
            motifIndex: -1,
            length: seq.sequence.length - previousEnd
          });
        }
      } else if (ranges.length === 0) {
        // If no ranges but sequence exists (even if empty), ensure sample appears
        stackData.push({
          sample: seq.name,
          start: 0,
          end: 0,
          motif: 'Interruption',
          motifIndex: -1,
          length: 0
        });
      }
    }
    
    // Process heatmap data - optimized with pre-allocated arrays
    const heatData: Array<{ sample: string; motif: string; count: number }> = [];
    const heatDataMap = new Map<string, number>(); // Use Map for faster lookups
    
    for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
      const seq = sequences[seqIdx];
      const motifIds = motifIdsList[seqIdx] || [];
      const motifCounts: { [motif: string]: number } = {};
      
      // Optimize: use for loop instead of forEach
      for (let i = 0; i < motifIds.length; i++) {
        const motifId = motifIds[i];
        if (motifId && motifId !== '.' && motifId !== '') {
          const motifIdx = parseInt(motifId);
          if (!isNaN(motifIdx) && motifIdx >= 0 && motifIdx < record.motifs.length) {
            const motifName = record.motifs[motifIdx];
            motifCounts[motifName] = (motifCounts[motifName] || 0) + 1;
          }
        }
      }
      
      // Use Object.entries but push directly
      for (const [motif, count] of Object.entries(motifCounts)) {
        heatData.push({ sample: seq.name, motif, count });
      }
    }
  
    // Get unique motifs and sort - optimized
    const motifTotalsMap = new Map<string, number>();
    for (const d of heatData) {
      motifTotalsMap.set(d.motif, (motifTotalsMap.get(d.motif) || 0) + d.count);
    }
    
    const motifTotals = Array.from(motifTotalsMap.entries())
      .map(([motif, total]) => ({ motif, total }))
      .sort((a, b) => b.total - a.total);
    const sortedMots = motifTotals.map(m => m.motif);
    
    // Color mapping - optimized: use Set directly and single pass
    const uniqueStackMotifsSet = new Set<string>();
    for (const d of stackData) {
      if (d.motif !== 'Interruption') {
        uniqueStackMotifsSet.add(d.motif);
      }
    }
    const uniqueStackMotifs = Array.from(uniqueStackMotifsSet);
    const colorMap: { [key: string]: string } = {};
    for (let idx = 0; idx < uniqueStackMotifs.length; idx++) {
      const motif = uniqueStackMotifs[idx];
      const motifIdx = record.motifs.indexOf(motif);
      if (motifIdx >= 0) {
        colorMap[motif] = COLOR_PALETTE[motifIdx % COLOR_PALETTE.length];
      } else {
        colorMap[motif] = COLOR_PALETTE[idx % COLOR_PALETTE.length];
      }
    }
    colorMap['Interruption'] = '#ef4444'; // Red color for interruptions
    
    // Optimize: calculate max in single pass instead of map + spread
    let maxLen = 1;
    let maxCnt = 1;
    for (const d of stackData) {
      if (d.end > maxLen) maxLen = d.end;
    }
    for (const d of heatData) {
      if (d.count > maxCnt) maxCnt = d.count;
    }
    
    // Pre-compute sample data for faster rendering (avoid filtering in render loop)
    // Optimize: build map directly instead of filtering
    const sampleDataMap: { [sample: string]: Array<{
      start: number;
      end: number;
      motif: string;
      motifIndex: number;
      length: number;
    }> } = {};
    
    // Initialize all samples with empty arrays
    for (const sample of sortedSamples) {
      sampleDataMap[sample] = [];
    }
    
    // Build map in one pass
    for (const item of stackData) {
      if (sampleDataMap[item.sample]) {
        sampleDataMap[item.sample].push(item);
      }
    }
    
    // Sort each sample's data once
    for (const sample of sortedSamples) {
      sampleDataMap[sample].sort((a, b) => a.start - b.start);
    }
    
    const result = {
      motifData: stackData,
      heatmapData: heatData,
      motifColorMap: colorMap,
      maxLength: maxLen,
      maxCount: maxCnt,
      sortedMotifs: sortedMots,
      sampleDataMap: sampleDataMap // Pre-computed per-sample data
    };
    
    // Cache the result
    if (plotCacheKey) {
      combinedPlotCache.set(plotCacheKey, result);
      
      // Limit cache size to prevent memory issues (keep last 50 entries)
      if (combinedPlotCache.size > 50) {
        const firstKey = combinedPlotCache.keys().next().value;
        if (firstKey) {
          combinedPlotCache.delete(firstKey);
        }
      }
    }
    
    return result;
  }, [sequences, spanList, motifIdsList, record.motifs, plotCacheKey, sortedSamples]);
  
  // Destructure plotData - use empty defaults if plotData is null/undefined
  const { motifData = [], heatmapData = [], motifColorMap = {}, maxLength = 1, maxCount = 1, sortedMotifs = [], sampleDataMap = {} } = plotData || {};
  
  // Create heatmap data map - MUST be called before any early returns
  const heatmapDataMap = useMemo(() => {
    const map: { [key: string]: number } = {};
    heatmapData.forEach(d => {
      map[`${d.sample}_${d.motif}`] = d.count;
    });
    return map;
  }, [heatmapData]);
  
  // Color function for heatmap - MUST be called before any early returns
  const getHeatmapColor = useMemo(() => {
    return (count: number) => {
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
  }, [maxCount]);
  
  // Memoize unique motifs for legend - MUST be called before any early returns
  const uniqueMotifsForLegend = useMemo(() => {
    const uniqueMotifs = new Set<string>();
    for (const d of motifData) {
      if (d.motif !== 'Interruption') {
        uniqueMotifs.add(d.motif);
      }
    }
    return Array.from(uniqueMotifs);
  }, [motifData]);
  
  // Early return check AFTER all hooks
  if (!plotData || !sequences || sequences.length === 0 || !sortedSamples || sortedSamples.length === 0) {
    return <div className="info-message">No data to display</div>;
  }
  
  return (
    <div className="combined-stack-heatmap">
      {/* Motif Legend at the top */}
      <div className="combined-legend">
        <div className="legend-section">
          <div className="legend-title">Motifs:</div>
          <div className="legend-items-scrollable">
            {uniqueMotifsForLegend.map(motif => (
              <div key={motif} className="legend-item-compact">
                <span className="legend-color-compact" style={{ backgroundColor: motifColorMap[motif] }} />
                <span className="legend-text-compact" title={motif}>
                  {motif.length > 12 ? motif.substring(0, 12) + '...' : motif}
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
        {sortedSamples.map(sample => {
          // Use pre-computed sample data (instant lookup, no filtering/sorting)
          const sampleStackData = sampleDataMap[sample] || [];
          // Optimize: calculate max in single pass instead of map + spread
          let totalLength = 0;
          for (const d of sampleStackData) {
            if (d.end > totalLength) totalLength = d.end;
          }
          // Ensure minimum length of 1 to avoid division by zero
          if (totalLength === 0) totalLength = 1;
          
          return (
            <div 
              key={sample} 
              className={`combined-row ${sample === 'Ref' || sample.startsWith('Allel') ? 'current-sample-row' : ''}`}
            >
              <div className="combined-sample-label">{sample}</div>
              <div className="combined-heatmap-cells-wrapper">
                <div className="combined-heatmap-cells">
                  {sortedMotifs.map(motif => {
                    const count = heatmapDataMap[`${sample}_${motif}`] || 0;
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
              <div className="combined-stack-bars-container" style={{ position: 'relative' }}>
                {/* Pathogenic threshold line */}
                {pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0 && (() => {
                  const thresholdNum = typeof pathogenicThreshold === 'number' ? pathogenicThreshold : parseFloat(String(pathogenicThreshold));
                  if (isNaN(thresholdNum) || thresholdNum <= 0) return null;
                  
                  const avgMotifLength = record.motifs.length > 0 ? record.motifs[0].length : 1;
                  const pathogenicThresholdLength = thresholdNum * avgMotifLength;
                  const isOutsideRange = pathogenicThresholdLength > maxLength;
                  const thresholdPercent = isOutsideRange 
                    ? 100 // Show at right edge if outside range
                    : (pathogenicThresholdLength / maxLength) * 100;
                  const exceedsThreshold = totalLength > pathogenicThresholdLength;
                  
                  return (
                    <>
                      <div
                        style={{
                          position: 'absolute',
                          left: `${thresholdPercent}%`,
                          top: 0,
                          bottom: 0,
                          width: '2px',
                          backgroundColor: '#EF4444',
                          borderLeft: isOutsideRange ? '2px solid #EF4444' : '2px dashed #EF4444',
                          zIndex: 10,
                          pointerEvents: 'none',
                          opacity: isOutsideRange ? 0.5 : 1
                        }}
                        title={`Pathogenic Threshold: ${thresholdNum} copies (${pathogenicThresholdLength}bp)${isOutsideRange ? ' (beyond visible range)' : ''}`}
                      />
                      {exceedsThreshold && (
                        <div
                          style={{
                            position: 'absolute',
                            left: `${Math.max((totalLength / maxLength) * 100, thresholdPercent)}%`,
                            top: '50%',
                            transform: 'translateY(-50%)',
                            marginLeft: '8px',
                            color: '#DC2626',
                            fontWeight: 'bold',
                            fontSize: '12px',
                            whiteSpace: 'nowrap',
                            zIndex: 11
                          }}
                          title="Exceeds pathogenic threshold"
                        >
                          🚨 PATHOGENIC
                        </div>
                      )}
                    </>
                  );
                })()}
                <div className="combined-stack-bars" style={{ width: `${Math.max((totalLength / maxLength) * 100, 0)}%` }}>
                  {sampleStackData.length > 0 ? sampleStackData.map((segment, idx) => {
                    const widthPercent = totalLength > 0 ? (segment.length / totalLength) * 100 : 0;
                    const color = motifColorMap[segment.motif] || '#9ca3af';
                    return (
                      <div
                        key={idx}
                        className="stack-plot-segment"
                        style={{
                          width: `${widthPercent}%`,
                          backgroundColor: color,
                          borderRight: idx < sampleStackData.length - 1 ? '1px solid rgba(0,0,0,0.1)' : 'none'
                        }}
                        title={`${segment.motif}: ${segment.start}-${segment.end} (${segment.length}bp)`}
                      />
                    );
                  }) : (
                    <div 
                      className="stack-plot-segment"
                      style={{
                        width: '100%',
                        backgroundColor: '#f3f4f6',
                        opacity: 0.5
                      }}
                      title="No data"
                    />
                  )}
                </div>
              </div>
              <div className="combined-length-label">{totalLength > 0 ? `${totalLength}bp` : '0bp'}</div>
            </div>
          );
        })}
      </div>
      
      {/* Motif headers at the bottom - Scrollable */}
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
});
CombinedStackHeatmap.displayName = 'CombinedStackHeatmap';

// Heatmap Component
interface HeatmapProps {
  record: Record;
  sequences: Array<{ name: string; sequence: string }>;
  spanList: string[];
  motifIdsList: string[][];
  sortedSamples?: string[];
}

const Heatmap: React.FC<HeatmapProps> = memo(({ record, sequences, spanList, motifIdsList, sortedSamples: providedSortedSamples }) => {
  // Create heatmap data: count of each motif per sample
  const heatmapData = useMemo(() => {
    const data: Array<{ sample: string; motif: string; count: number }> = [];
    
    sequences.forEach((seq, seqIdx) => {
      const motifIds = motifIdsList[seqIdx] || [];
      const motifCounts: { [motif: string]: number } = {};
      
      // Count each motif type
      motifIds.forEach(motifId => {
        if (motifId && motifId !== '.' && motifId !== '') {
          const motifIdx = parseInt(motifId);
          if (!isNaN(motifIdx) && motifIdx >= 0 && motifIdx < record.motifs.length) {
            const motifName = record.motifs[motifIdx];
            motifCounts[motifName] = (motifCounts[motifName] || 0) + 1;
          }
        }
      });
      
      // Add to data array
      Object.entries(motifCounts).forEach(([motif, count]) => {
        data.push({ sample: seq.name, motif, count });
      });
    });
    
    return data;
  }, [sequences, motifIdsList, record.motifs]);
  
  // Get unique motifs
  const uniqueMotifs = Array.from(new Set(heatmapData.map(d => d.motif)));
  
  // Use provided sorted samples if available, otherwise calculate
  const sortedSamples = useMemo(() => {
    if (providedSortedSamples && providedSortedSamples.length > 0) {
      return providedSortedSamples;
    }
    // Fallback: sort by total count
    const uniqueSamples = Array.from(new Set(heatmapData.map(d => d.sample)));
    const sampleTotals = uniqueSamples.map(sample => {
      const total = heatmapData.filter(d => d.sample === sample).reduce((sum, d) => sum + d.count, 0);
      return { sample, total };
    }).sort((a, b) => b.total - a.total);
    return sampleTotals.map(s => s.sample);
  }, [providedSortedSamples, heatmapData]);
  
  // Sort motifs by total count (descending)
  const motifTotals = uniqueMotifs.map(motif => {
    const total = heatmapData.filter(d => d.motif === motif).reduce((sum, d) => sum + d.count, 0);
    return { motif, total };
  }).sort((a, b) => b.total - a.total);
  
  const sortedMotifs = motifTotals.map(m => m.motif);
  
  // Get max count for color scaling
  const maxCount = Math.max(...heatmapData.map(d => d.count), 1);
  
  // Create color scale function (red-purple gradient)
  const getColor = (count: number) => {
    const intensity = count / maxCount;
    // Red-purple gradient: from light red to dark purple
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
  
  // Create a map for quick lookup
  const dataMap: { [key: string]: number } = {};
  heatmapData.forEach(d => {
    dataMap[`${d.sample}_${d.motif}`] = d.count;
  });
  
  return (
    <div className="heatmap-container">
      <div className="heatmap-chart">
        <div className="heatmap-header">
          {sortedMotifs.map(motif => (
            <div key={motif} className="heatmap-motif-header">{motif}</div>
          ))}
        </div>
        <div className="heatmap-body">
          {sortedSamples.map(sample => (
            <div key={sample} className="heatmap-row">
              {/* No sample label here - it's shared with stack plot */}
              <div className="heatmap-cells">
                {sortedMotifs.map(motif => {
                  const count = dataMap[`${sample}_${motif}`] || 0;
                  const color = getColor(count);
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
                      {count > 0 && <span className="heatmap-cell-value">{count}</span>}
                    </div>
                  );
                })}
              </div>
            </div>
          ))}
        </div>
      </div>
      
      {/* Color scale legend */}
      <div className="heatmap-legend">
        <div className="legend-label">Count:</div>
        <div className="legend-gradient">
          <span>0</span>
          <div className="gradient-bar" />
          <span>{maxCount}</span>
        </div>
      </div>
    </div>
  );
});
Heatmap.displayName = 'Heatmap';

// Bar Plot Component
interface BarPlotProps {
  sequences: Array<{ name: string; sequence: string }>;
  motifIdsList: string[][];
  pathogenicThreshold?: number; // Pathogenic threshold in copy number
}

const BarPlot: React.FC<BarPlotProps> = memo(({ sequences, motifIdsList, pathogenicThreshold }) => {
  // Count actual motifs per sample (excluding interruptions)
  const motifCounts = useMemo(() => {
    return sequences.map((seq, idx) => {
      const motifIds = motifIdsList[idx] || [];
      const validMotifs = motifIds.filter(id => id && id !== '.' && id !== '');
      return {
        sample: seq.name,
        count: validMotifs.length
      };
    }).sort((a, b) => b.count - a.count);
  }, [sequences, motifIdsList]);
  
  const maxCount = Math.max(...motifCounts.map(m => m.count), 1);
  
  // Extend maxCount to include threshold if it's above current max
  const effectiveMaxCount = useMemo(() => {
    if (pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0) {
      const thresholdNum = typeof pathogenicThreshold === 'number' ? pathogenicThreshold : parseFloat(String(pathogenicThreshold));
      if (!isNaN(thresholdNum) && thresholdNum > maxCount) {
        // Round up to next nice number (e.g., if threshold is 20, use 20 or round to 25)
        return Math.max(thresholdNum, Math.ceil(thresholdNum / 5) * 5);
      }
    }
    return maxCount;
  }, [maxCount, pathogenicThreshold]);
  
  // Calculate statistics
  const stats = useMemo(() => {
    const total = motifCounts.reduce((sum, item) => sum + item.count, 0);
    const avg = total / motifCounts.length;
    const sorted = [...motifCounts].sort((a, b) => a.count - b.count);
    const median = sorted.length > 0 
      ? sorted[Math.floor(sorted.length / 2)].count 
      : 0;
    return { total, avg: Math.round(avg * 10) / 10, median };
  }, [motifCounts]);
  
  // Generate y-axis ticks (0 to effectiveMaxCount)
  const yAxisTicks = useMemo(() => {
    const numTicks = 6;
    const ticks: number[] = [];
    for (let i = 0; i <= numTicks; i++) {
      const value = Math.round((effectiveMaxCount / numTicks) * i);
      ticks.push(value);
    }
    return ticks.reverse(); // Highest to lowest (top to bottom)
  }, [effectiveMaxCount]);
  
  const CHART_HEIGHT = 350;
  const BAR_WIDTH = 60;
  const BAR_GAP = 20;
  const Y_AXIS_WIDTH = 60;
  
  return (
    <div className="bar-plot-v2">
      {/* Statistics Cards */}
      <div className="bar-plot-stats-v2">
        <div className="stat-card-v2">
          <div className="stat-info-v2">
            <div className="stat-label-v2">Total Motifs</div>
            <div className="stat-value-v2">{stats.total}</div>
          </div>
        </div>
        <div className="stat-card-v2">
          <div className="stat-info-v2">
            <div className="stat-label-v2">Average</div>
            <div className="stat-value-v2">{stats.avg}</div>
          </div>
        </div>
        <div className="stat-card-v2">
          <div className="stat-info-v2">
            <div className="stat-label-v2">Median</div>
            <div className="stat-value-v2">{stats.median}</div>
          </div>
        </div>
        <div className="stat-card-v2">
          <div className="stat-info-v2">
            <div className="stat-label-v2">Maximum</div>
            <div className="stat-value-v2">{maxCount}</div>
          </div>
        </div>
      </div>
      
      {/* Chart with unified coordinate system */}
      <div className="bar-plot-chart-new">
        <div className="chart-unified-container" style={{ height: `${CHART_HEIGHT}px` }}>
          {/* Y-axis labels - positioned with exact pixel calculations */}
          <div className="y-axis-new" style={{ width: `${Y_AXIS_WIDTH}px`, height: `${CHART_HEIGHT}px` }}>
            <div className="y-axis-label-new">Count</div>
            {yAxisTicks.map((tick, i) => {
              const yPosition = (i / (yAxisTicks.length - 1)) * CHART_HEIGHT;
              return (
                <div 
                  key={i} 
                  className="y-tick-new"
                  style={{ top: `${yPosition}px` }}
                >
                  <span className="tick-line-new"></span>
                  <span className="tick-label-new">{tick}</span>
                </div>
              );
            })}
          </div>
          
          {/* Chart area with grid and bars */}
          <div className="chart-area-new" style={{ flex: 1, height: `${CHART_HEIGHT}px`, position: 'relative' }}>
            {/* Grid lines - using same pixel calculations */}
            {yAxisTicks.map((tick, i) => {
              const yPosition = (i / (yAxisTicks.length - 1)) * CHART_HEIGHT;
              return (
                <div 
                  key={i} 
                  className="grid-line-new"
                  style={{ top: `${yPosition}px` }}
                ></div>
              );
            })}
            
            {/* Pathogenic threshold line */}
            {pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0 && (() => {
              const thresholdNum = typeof pathogenicThreshold === 'number' ? pathogenicThreshold : parseFloat(String(pathogenicThreshold));
              if (isNaN(thresholdNum) || thresholdNum <= 0) return null;
              
              // Calculate position based on Y-axis scale (0 to effectiveMaxCount)
              // The Y-axis represents counts from 0 (bottom) to effectiveMaxCount (top)
              // Position threshold at its actual value on the Y-axis
              const bottomPosition = (thresholdNum / effectiveMaxCount) * CHART_HEIGHT;
              
              return (
                <>
                  <div
                    style={{
                      position: 'absolute',
                      left: 0,
                      right: 0,
                      bottom: `${bottomPosition}px`,
                      height: '2px',
                      backgroundColor: '#EF4444',
                      borderTop: '2px dashed #EF4444',
                      zIndex: 10,
                      pointerEvents: 'none'
                    }}
                    title={`Pathogenic Threshold: ${thresholdNum} repeats`}
                  />
                  <div
                    style={{
                      position: 'absolute',
                      left: '8px',
                      bottom: `${bottomPosition - 20}px`,
                      color: '#DC2626',
                      fontWeight: 'bold',
                      fontSize: '16px',
                      zIndex: 11,
                      pointerEvents: 'none',
                      backgroundColor: 'rgba(255, 255, 255, 0.9)',
                      padding: '2px 6px',
                      borderRadius: '4px'
                    }}
                  >
                    🚨 Pathogenic Threshold
                  </div>
                </>
              );
            })()}
            
            {/* Bars */}
            <div className="bars-container-new" style={{ paddingLeft: '0.5rem' }}>
              {motifCounts.map((item, index) => {
                const barHeight = effectiveMaxCount > 0 ? (item.count / effectiveMaxCount) * CHART_HEIGHT : 0;
                const color = COLOR_PALETTE[index % COLOR_PALETTE.length];
                const leftPosition = index * (BAR_WIDTH + BAR_GAP);
                const exceedsThreshold = pathogenicThreshold && pathogenicThreshold > 0 && item.count >= pathogenicThreshold;
                
                return (
                  <div
                    key={item.sample}
                    className="bar-new"
                    style={{
                      position: 'absolute',
                      left: `${leftPosition}px`,
                      bottom: '0',
                      width: `${BAR_WIDTH}px`,
                      height: `${barHeight}px`,
                      backgroundColor: color,
                      border: exceedsThreshold ? '3px solid #DC2626' : 'none',
                      boxShadow: exceedsThreshold ? '0 0 8px rgba(220, 38, 38, 0.5)' : 'none'
                    }}
                    title={`${item.sample}: ${item.count} motifs${exceedsThreshold ? ' (PATHOGENIC)' : ''}`}
                  >
                    {barHeight > 30 && (
                      <span className="bar-count-new">{item.count}</span>
                    )}
                  </div>
                );
              })}
            </div>
          </div>
        </div>
        
        {/* X-axis labels */}
        <div className="x-axis-labels-new">
          <div className="x-axis-labels-container" style={{ paddingLeft: `${Y_AXIS_WIDTH + 0.5}rem` }}>
            {motifCounts.map((item, index) => {
              const leftPosition = index * (BAR_WIDTH + BAR_GAP);
              return (
                <div
                  key={item.sample}
                  className="x-label-new"
                  style={{ 
                    left: `${leftPosition}px`,
                    width: `${BAR_WIDTH}px`
                  }}
                >
                  {item.sample}
                </div>
              );
            })}
          </div>
        </div>
      </div>
    </div>
  );
});
BarPlot.displayName = 'BarPlot';

// Cluster Plot Component (simplified - shows sample stats for clustering)
interface ClusterPlotProps {
  sequences: Array<{ name: string; sequence: string }>;
  motifIdsList: string[][];
}

const ClusterPlot: React.FC<ClusterPlotProps> = memo(({ sequences, motifIdsList }) => {
  // Calculate copy number and length for each sample
  const sampleStats = useMemo(() => {
    return sequences.map((seq, idx) => {
      const motifIds = motifIdsList[idx] || [];
      const validMotifs = motifIds.filter(id => id && id !== '.' && id !== '');
      const copyNumber = validMotifs.length;
      const length = seq.sequence.length;
      
      return {
        sample: seq.name,
        copyNumber,
        length
      };
    });
  }, [sequences, motifIdsList]);
  
  const maxCopyNumber = Math.max(...sampleStats.map(s => s.copyNumber), 1);
  const maxLength = Math.max(...sampleStats.map(s => s.length), 1);
  
  // Simple clustering visualization: scatter plot of copy number vs length
  return (
    <div className="cluster-plot-container">
      <div className="cluster-plot-scatter">
        <div className="cluster-plot-axis-labels">
          <div className="y-axis-label">Copy Number</div>
          <div className="cluster-plot-area">
            {sampleStats.map((stat, idx) => {
              const xPercent = (stat.length / maxLength) * 100;
              const yPercent = 100 - (stat.copyNumber / maxCopyNumber) * 100; // Invert Y for visual
              const color = COLOR_PALETTE[idx % COLOR_PALETTE.length];
              const isCurrent = stat.sample === 'Current Sample' || stat.sample.startsWith('Allel') || stat.sample === 'Ref';
              
              return (
                <div
                  key={stat.sample}
                  className={`cluster-plot-point ${isCurrent ? 'current-sample' : ''}`}
                  style={{
                    left: `${xPercent}%`,
                    bottom: `${yPercent}%`,
                    backgroundColor: color,
                    borderColor: isCurrent ? '#dc2626' : color,
                    borderWidth: isCurrent ? '3px' : '2px'
                  }}
                  title={`${stat.sample}: Copy Number=${stat.copyNumber}, Length=${stat.length}bp`}
                >
                  <span className="cluster-point-label">{stat.sample}</span>
                </div>
              );
            })}
            
            {/* Grid lines */}
            <div className="cluster-grid-lines">
              {[0, 25, 50, 75, 100].map(percent => (
                <React.Fragment key={percent}>
                  <div className="grid-line-vertical" style={{ left: `${percent}%` }} />
                  <div className="grid-line-horizontal" style={{ bottom: `${percent}%` }} />
                </React.Fragment>
              ))}
            </div>
          </div>
          <div className="x-axis-label">Length (bp)</div>
        </div>
      </div>
      
      <div className="cluster-plot-legend">
        <div className="legend-note">
          <span className="legend-current-marker" /> = Current Sample / Reference
        </div>
      </div>
      
      <div className="cluster-plot-stats">
        <h4>Sample Statistics</h4>
        <div className="cluster-stats-grid">
          {sampleStats.map(stat => (
            <div key={stat.sample} className="cluster-stat-item">
              <div className="stat-sample">{stat.sample}</div>
              <div className="stat-values">
                <span>Copy #: {stat.copyNumber}</span>
                <span>Length: {stat.length}bp</span>
              </div>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
});
ClusterPlot.displayName = 'ClusterPlot';

export default memo(PopulationComparison);

