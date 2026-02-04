import React, { useState, useEffect, useMemo, memo, useCallback, useRef } from 'react';
import axios from 'axios';
import './PopulationComparison.css'; // Reuse the same styles
import ExportMenu from './ExportMenu';
import RegionInfoCard from './RegionInfoCard';
import { exportToCSV, exportToFASTA, generateFilename } from '../utils/exportUtils';

// Import visualization components from PopulationComparison
// We'll need to make these accessible or duplicate the logic
// For now, we'll create a simplified version that uses the same components

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

// Reuse color palette from PopulationComparison
const COLOR_PALETTE = [
  '#667eea', '#764ba2', '#f093fb', '#f5576c', '#4fd1c7', '#68d391',
  '#f6e05e', '#f6ad55', '#fc8181', '#7e9af9', '#c084fc', '#f472b6',
  '#60a5fa', '#34d399', '#fbbf24', '#fb7185', '#a78bfa', '#38bdf8'
];

// Cache for combined plot data (shared with PopulationComparison if needed)
const combinedPlotCache = new Map<string, {
  motifData: Array<{ sample: string; start: number; end: number; motif: string; motifIndex: number; length: number }>;
  heatmapData: Array<{ sample: string; motif: string; count: number }>;
  motifColorMap: { [key: string]: string };
  maxLength: number;
  maxCount: number;
  sortedMotifs: string[];
  sampleDataMap: { [sample: string]: Array<{ start: number; end: number; motif: string; motifIndex: number; length: number }> };
}>();

// Cache for parsed motif ranges
const parsedRangeCache = new Map<string, Array<[number, number]>>();

interface PopulationRecord {
  chr: string;
  pos: number;
  stop: number;
  motifs: string[];
  motif_ids_h?: string[];  // For assembly/single-haplotype records
  motif_ids_h1?: string[]; // For diploid records
  motif_ids_h2?: string[];  // For diploid records
  motif_ids_ref: string[];
  ref_CN: number;
  CN_H?: number;  // For assembly/single-haplotype records
  CN_H1?: string | null;  // For diploid records
  CN_H2?: string | null;  // For diploid records
  spans: string | string[];
  ref_allele: string;
  alt_allele?: string;  // For assembly/single-haplotype records
  alt_allele1?: string;  // For diploid records
  alt_allele2?: string;  // For diploid records
  gt: string;  // Can be "0", "1" (single-allele) or "0/0", "0/1", "1/1" (diploid)
  id: string;
  // Additional fields for haplotype-specific files
  original_sample_name?: string;
  base_sample_name?: string;
  haplotype?: string;
}

type CohortMode = 'cohort-read' | 'cohort-assembly';

interface CohortAnalysisProps {
  mode: CohortMode;
  region: string;
  publicVcfFolder?: string;
  record?: any; // Optional record for reference/allele comparison
}

const CohortAnalysis: React.FC<CohortAnalysisProps> = ({ mode, region, publicVcfFolder, record }) => {
  const [populationRecords, setPopulationRecords] = useState<{ [key: string]: PopulationRecord }>({});
  const [sampleIds, setSampleIds] = useState<Array<{ sample_name: string; file_path: string }>>([]);
  const [loadedSamples, setLoadedSamples] = useState<Set<string>>(new Set());
  const [loading, setLoading] = useState(false);
  const [loadingIds, setLoadingIds] = useState(false);
  const [loadingRecords, setLoadingRecords] = useState(false);
  const [error, setError] = useState<string | null>(null);
  const [firstSampleRecord, setFirstSampleRecord] = useState<PopulationRecord | null>(null);
  const [selectedSample, setSelectedSample] = useState<string | null>(null);
  const [searchQuery, setSearchQuery] = useState('');
  const [selectedGenotypes, setSelectedGenotypes] = useState<string[]>([]);
  const [sortBy, setSortBy] = useState<'length' | 'motifCount' | 'name' | 'genotype'>('length');
  const [sortDirection, setSortDirection] = useState<'asc' | 'desc'>('desc');
  const [pathogenicInfo, setPathogenicInfo] = useState<{ 
    pathogenic_threshold?: number;
    gene?: string;
    disease?: string;
    inheritance?: string;
  } | null>(null);
  
  // Batch size for loading records - increased for better throughput
  const BATCH_SIZE = 50;
  
  // Collapsible sections state
  const [expandedSections, setExpandedSections] = useState({
    frequency: true,
    genotype: true,
    stackHeatmap: true,
    barPlot: true,
    cluster: true,
    samples: true
  });
  // Quick filter for sample table: 'all' | 'loaded' | 'not-loaded'
  const [statusFilter, setStatusFilter] = useState<'all' | 'loaded' | 'not-loaded'>('all');

  // Virtual scrolling state
  const [virtualScrollEnabled, setVirtualScrollEnabled] = useState(false);
  const [visibleRange, setVisibleRange] = useState({ start: 0, end: 50 });
  const VIRTUAL_SCROLL_THRESHOLD = 50;

  // Ref to prevent concurrent loadRecordsBatch calls (set synchronously before async work)
  const loadingRecordsRef = useRef(false);

  // Create sequences data from population records
  const sequencesData = useMemo(() => {
    if (!publicVcfFolder || Object.keys(populationRecords).length === 0) {
      return { sequences: [], spanList: [], motifIdsList: [], sortedSamples: [] };
    }

    const sequences: Array<{ name: string; sequence: string }> = [];
    const spanList: string[] = [];
    const motifIdsList: string[][] = [];

    // Extract reference from first available record (all records have the same ref_allele for the same region)
    // Priority: 1) record prop, 2) firstSampleRecord, 3) first loaded population record
    let refData: { ref_allele: string; spans: string | string[]; motif_ids_ref: string[] } | null = null;
    
    if (record && record.ref_allele) {
      // Use record prop if available
      refData = {
        ref_allele: record.ref_allele,
        spans: record.spans || [],
        motif_ids_ref: record.motif_ids_ref || []
      };
    } else if (firstSampleRecord && firstSampleRecord.ref_allele) {
      // Use first sample record from backend
      refData = {
        ref_allele: firstSampleRecord.ref_allele,
        spans: firstSampleRecord.spans || '',
        motif_ids_ref: firstSampleRecord.motif_ids_ref || []
      };
    } else {
      // Extract from first loaded population record
      const firstRecord = Object.values(populationRecords)[0];
      if (firstRecord && firstRecord.ref_allele) {
        refData = {
          ref_allele: firstRecord.ref_allele,
          spans: firstRecord.spans || '',
          motif_ids_ref: firstRecord.motif_ids_ref || []
        };
      }
    }
    
    // Add reference to sequences if available
    if (refData && refData.ref_allele) {
      sequences.push({ name: 'Ref', sequence: refData.ref_allele });
      const spansValue = Array.isArray(refData.spans) ? refData.spans[0] || '' : (refData.spans || '');
      spanList.push(spansValue);
      motifIdsList.push(refData.motif_ids_ref || []);
    }

    // Process population samples - now using sample names from VCF headers
    // Samples are no longer grouped by h1/h2 suffix, just use the sample name directly
    for (const [sampleName, popRecord] of Object.entries(populationRecords)) {
      // Check if this record has alt_allele (assembly format) or alt_allele1/alt_allele2 (regular format)
      const recordAny = popRecord as any;
      
      // #region agent log
      fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:118','message':'processing sample record',data:{sampleName,has_alt_allele:!!popRecord.alt_allele,has_alt_allele1:recordAny.alt_allele1!==undefined,has_alt_allele2:recordAny.alt_allele2!==undefined,alt_allele:popRecord.alt_allele?.substring(0,50)||'',alt_allele1:recordAny.alt_allele1?.substring(0,50)||'',alt_allele2:recordAny.alt_allele2?.substring(0,50)||'',alt_allele1_empty:!recordAny.alt_allele1||recordAny.alt_allele1==='',alt_allele2_empty:!recordAny.alt_allele2||recordAny.alt_allele2===''},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      
      if (popRecord.alt_allele && popRecord.alt_allele !== '') {
        // Assembly format - single haplotype
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:123','message':'using assembly format',data:{sampleName},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
        // #endregion
        sequences.push({ name: sampleName, sequence: popRecord.alt_allele });
        // Handle spans as either string or string array
        const spansValue = Array.isArray(popRecord.spans) ? popRecord.spans[0] || '' : (popRecord.spans || '');
        spanList.push(spansValue);
        motifIdsList.push(popRecord.motif_ids_h || []);
      } else if (recordAny.alt_allele1 !== undefined) {
        // Regular format - may have multiple haplotypes
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:129','message':'using diploid format',data:{sampleName,alt_allele1_exists:!!recordAny.alt_allele1,alt_allele1_nonempty:recordAny.alt_allele1&&recordAny.alt_allele1!=='',alt_allele2_exists:!!recordAny.alt_allele2,alt_allele2_nonempty:recordAny.alt_allele2&&recordAny.alt_allele2!==''},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
        // #endregion
        if (recordAny.alt_allele1 && recordAny.alt_allele1 !== '') {
          // #region agent log
          fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:131','message':'adding h1 sequence',data:{sampleName,sequenceLength:recordAny.alt_allele1.length,sequencePreview:recordAny.alt_allele1.substring(0,50)},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
          // #endregion
          sequences.push({ name: `${sampleName}_h1`, sequence: recordAny.alt_allele1 });
          spanList.push(Array.isArray(recordAny.spans) ? recordAny.spans[1] || '' : (recordAny.spans || ''));
          motifIdsList.push(recordAny.motif_ids_h1 || []);
        }
        if (recordAny.alt_allele2 && recordAny.alt_allele2 !== '') {
          // #region agent log
          fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:136','message':'adding h2 sequence',data:{sampleName,sequenceLength:recordAny.alt_allele2.length,sequencePreview:recordAny.alt_allele2.substring(0,50)},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
          // #endregion
          sequences.push({ name: `${sampleName}_h2`, sequence: recordAny.alt_allele2 });
          spanList.push(Array.isArray(recordAny.spans) ? recordAny.spans[2] || '' : (recordAny.spans || ''));
          motifIdsList.push(recordAny.motif_ids_h2 || []);
        }
      } else {
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:141','message':'no allele format matched',data:{sampleName},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
        // #endregion
      }
    }

    // Sort samples by actual sequence length (not span-based length)
    const sampleLengths = sequences.map((seq) => {
      const totalLength = seq.sequence ? seq.sequence.length : 0;
      return { sample: seq.name, totalLength };
    }).sort((a, b) => b.totalLength - a.totalLength);
    
    const sortedSamples = sampleLengths.map(s => s.sample);

    console.log('CohortAnalysis: sequencesData created', {
      sequenceCount: sequences.length,
      sampleNames: sequences.map(s => s.name).slice(0, 5),
      sortedSamples: sortedSamples.slice(0, 5)
    });

    return { sequences, spanList, motifIdsList, sortedSamples };
  }, [populationRecords, publicVcfFolder, record, firstSampleRecord]);

  // Genotype comparison data
  const genotypeComparisonData = useMemo(() => {
    if (Object.keys(populationRecords).length === 0) return {};
    
    const genotypes: { [sample: string]: string } = {};
    
    // Add reference genotype - use '0/0' as default (homozygous reference)
    // or extract from record if available
    if (record && record.gt) {
      genotypes['Ref'] = record.gt;
    } else if (firstSampleRecord && firstSampleRecord.gt) {
      genotypes['Ref'] = firstSampleRecord.gt;
    } else {
      // Default to homozygous reference genotype
      genotypes['Ref'] = '0/0';
    }

    // Process samples - check if they have regular format (alt_allele1/alt_allele2) or assembly format (alt_allele)
    const sampleGroups: { [base: string]: { h1?: any; h2?: any } } = {};
    
    for (const [sampleName, popRecord] of Object.entries(populationRecords)) {
      const recordAny = popRecord as any;
      
      // Check if it's regular format (has alt_allele1) or assembly format (has alt_allele)
      if (recordAny.alt_allele1 !== undefined) {
        // Regular format - may have h1 and h2 haplotypes
        // Group by sample name and process separately
        if (!sampleGroups[sampleName]) sampleGroups[sampleName] = {};
        if (recordAny.alt_allele1) {
          sampleGroups[sampleName].h1 = {
            // Handle both single-allele (0, 1) and diploid (0/0, 0/1) genotypes
            gt: recordAny.gt?.includes('/') ? recordAny.gt.split('/')[0] : (recordAny.gt || '0'),
            motif_ids_h: recordAny.motif_ids_h1 || []
          };
        }
        if (recordAny.alt_allele2) {
          sampleGroups[sampleName].h2 = {
            // Handle both single-allele (0, 1) and diploid (0/0, 0/1) genotypes
            gt: recordAny.gt?.includes('/') ? recordAny.gt.split('/')[1] : (recordAny.gt || '0'),
            motif_ids_h: recordAny.motif_ids_h2 || []
          };
        }
        // Use the full genotype string if available
        if (recordAny.gt) {
          genotypes[sampleName] = recordAny.gt;
        }
      } else {
        // Assembly format - single haplotype, use sample name directly
        genotypes[sampleName] = popRecord.gt || '0';
      }
    }
    
    // For samples with h1/h2 haplotypes, compute diploid genotype if needed
    for (const [baseName, group] of Object.entries(sampleGroups)) {
      if (group.h1 && group.h2) {
        const gt_h1 = group.h1.gt || '0';
        const gt_h2 = group.h2.gt || '0';
        const ids_h1 = group.h1.motif_ids_h || [];
        const ids_h2 = group.h2.motif_ids_h || [];
        // Only compute if we don't already have a combined genotype
        if (!genotypes[baseName]) {
          genotypes[baseName] = computeDiploidGenotype(gt_h1, gt_h2, ids_h1, ids_h2);
        }
      } else if (group.h1 && !genotypes[baseName]) {
        genotypes[baseName] = group.h1.gt || '0';
      } else if (group.h2 && !genotypes[baseName]) {
        genotypes[baseName] = group.h2.gt || '0';
      }
    }
    
    return genotypes;
  }, [populationRecords, record, firstSampleRecord]);

  // Fetch pathogenic info for the region
  useEffect(() => {
    const fetchPathogenicInfo = async () => {
      if (!region || !firstSampleRecord) {
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
  }, [region, firstSampleRecord]);
  
  // Get all sample names (from IDs, including not yet loaded)
  const allSampleNames = useMemo(() => {
    const names = sampleIds.map(id => id.sample_name);
    // Add reference if record is provided
    if (record) {
      names.unshift('Ref');
    }
    return names;
  }, [sampleIds, record]);

  // Filter and sort samples (works with both loaded and unloaded samples)
  const filteredSamples = useMemo(() => {
    let samples = allSampleNames;
    
    if (searchQuery) {
      samples = samples.filter(s => s.toLowerCase().includes(searchQuery.toLowerCase()));
    }
    
    if (selectedGenotypes.length > 0) {
      samples = samples.filter(sample => {
        // Only filter by genotype if the sample is loaded
        if (!loadedSamples.has(sample) && sample !== 'Ref') return true;
        const genotype = genotypeComparisonData[sample] || '';
        return selectedGenotypes.includes(String(genotype).trim());
      });
    }
    
    if (statusFilter === 'loaded') {
      samples = samples.filter(s => loadedSamples.has(s) || s === 'Ref');
    } else if (statusFilter === 'not-loaded') {
      samples = samples.filter(s => !loadedSamples.has(s) && s !== 'Ref');
    }
    
    // Sort - prioritize loaded samples, then sort by criteria
    const loaded = samples.filter(s => loadedSamples.has(s) || s === 'Ref');
    const unloaded = samples.filter(s => !loadedSamples.has(s) && s !== 'Ref');
    
    const sortedLoaded = [...loaded].sort((a, b) => {
      let comparison = 0;
      switch(sortBy) {
        case 'length':
          // For length, find the maximum length across all haplotypes for this sample
          const sequencesA = sequencesData.sequences.filter(s => {
            const baseName = s.name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
            return baseName === a || s.name === a;
          });
          const sequencesB = sequencesData.sequences.filter(s => {
            const baseName = s.name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
            return baseName === b || s.name === b;
          });
          const maxLenA = Math.max(...sequencesA.map(s => s.sequence.length), 0);
          const maxLenB = Math.max(...sequencesB.map(s => s.sequence.length), 0);
          comparison = maxLenA - maxLenB;
          break;
        case 'motifCount':
          // For motif count, sum across all haplotypes for this sample
          const seqsA = sequencesData.sequences.filter(s => {
            const baseName = s.name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
            return baseName === a || s.name === a;
          });
          const seqsB = sequencesData.sequences.filter(s => {
            const baseName = s.name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
            return baseName === b || s.name === b;
          });
          let totalMotifsA = 0;
          let totalMotifsB = 0;
          seqsA.forEach(seq => {
            const idx = sequencesData.sequences.findIndex(s => s.name === seq.name);
            if (idx >= 0) {
              const motifIds = sequencesData.motifIdsList[idx] || [];
              totalMotifsA += motifIds.filter(id => id && id !== '.' && id !== '').length;
            }
          });
          seqsB.forEach(seq => {
            const idx = sequencesData.sequences.findIndex(s => s.name === seq.name);
            if (idx >= 0) {
              const motifIds = sequencesData.motifIdsList[idx] || [];
              totalMotifsB += motifIds.filter(id => id && id !== '.' && id !== '').length;
            }
          });
          comparison = totalMotifsA - totalMotifsB;
          break;
        case 'genotype':
          const gtA = genotypeComparisonData[a] || '';
          const gtB = genotypeComparisonData[b] || '';
          comparison = String(gtA).localeCompare(String(gtB));
          break;
        case 'name':
          comparison = a.localeCompare(b);
          break;
      }
      return sortDirection === 'asc' ? comparison : -comparison;
    });
    
    return [...sortedLoaded, ...unloaded];
  }, [allSampleNames, sequencesData, searchQuery, selectedGenotypes, statusFilter, sortBy, sortDirection, genotypeComparisonData, loadedSamples]);

  // Virtual scrolling
  const visibleSamples = useMemo(() => {
    if (filteredSamples.length <= VIRTUAL_SCROLL_THRESHOLD) {
      return filteredSamples;
    }
    
    const refSamples = filteredSamples.filter(s => s === 'Ref' || s.startsWith('Allel'));
    const populationSamples = filteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel'));
    const visiblePopulationSamples = populationSamples.slice(visibleRange.start, visibleRange.end);
    
    return [...refSamples, ...visiblePopulationSamples];
  }, [filteredSamples, visibleRange]);

  useEffect(() => {
    setVirtualScrollEnabled(filteredSamples.length > VIRTUAL_SCROLL_THRESHOLD);
  }, [filteredSamples.length]);

  // Fetch population data
  // First, load just the sample IDs (fast)
  useEffect(() => {
    if (!publicVcfFolder || !region) {
      setSampleIds([]);
      setPopulationRecords({});
      setLoadedSamples(new Set());
      setLoadingIds(false);
      setFirstSampleRecord(null);
      setSelectedSample(null);
      setPathogenicInfo(null);
      return;
    }
    
    // Reset state when region changes
    setSampleIds([]);
    setPopulationRecords({});
    setLoadedSamples(new Set());
    setFirstSampleRecord(null);
    setSelectedSample(null);
    setPathogenicInfo(null);
    setError(null);

    const fetchSampleIds = async () => {
      setLoadingIds(true);
      setError(null);
      const startTime = Date.now();
      
      // #region agent log
      fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:313',message:'fetchSampleIds start',data:{region,publicVcfFolder,startTime},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
      // #endregion
      
      try {
        const url = `${API_BASE}/api/population/region/${encodeURIComponent(region)}/ids`;
        const apiStartTime = Date.now();
        
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:320',message:'fetchSampleIds API call start',data:{url},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
        
        const response = await axios.get(url, {
          params: { folder_path: publicVcfFolder, mode: mode }
        });
        
        const apiTime = Date.now() - apiStartTime;
        
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:328',message:'fetchSampleIds API response',data:{apiTime,status:response.status,sample_ids_count:response.data?.sample_ids?.length || 0,has_first_record:!!response.data?.first_record},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
        
        const ids = response.data.sample_ids || [];
        setSampleIds(ids);
        
        // If first sample's record is included, add it immediately
        if (response.data.first_record && response.data.first_sample) {
          const firstSampleName = response.data.first_sample;
          setFirstSampleRecord(response.data.first_record);
          setPopulationRecords(prev => ({
            ...prev,
            [firstSampleName]: response.data.first_record
          }));
          setLoadedSamples(prev => {
            const newSet = new Set(prev);
            newSet.add(firstSampleName);
            return newSet;
          });
        }
        
        // Load next batch of records (skip first one if already loaded)
        if (ids.length > 1) {
          const startIdx = response.data.first_record ? 1 : 0;
          const firstBatch = ids.slice(startIdx, startIdx + BATCH_SIZE)
            .map((id: { sample_name: string; file_path: string }) => id.sample_name);
          
          if (firstBatch.length > 0) {
            setTimeout(() => {
              loadRecordsBatch(firstBatch);
            }, 100);
          }
        }
        
        const totalTime = Date.now() - startTime;
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:360',message:'fetchSampleIds complete',data:{totalTime,apiTime},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
      } catch (err: any) {
        // #region agent log
        fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:364',message:'fetchSampleIds error',data:{error:err.message,status:err.response?.status},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
        console.error('Error loading sample IDs:', err);
        setError(err.response?.data?.detail || err.message || 'Failed to load sample IDs');
        setSampleIds([]);
      } finally {
        setLoadingIds(false);
      }
    };

    fetchSampleIds();
  }, [publicVcfFolder, region, mode]);

  // Function to load records for specific samples
  const loadRecordsBatch = useCallback(async (sampleNames: string[]) => {
    // #region agent log
    fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:343',message:'loadRecordsBatch entry',data:{sampleNamesCount:sampleNames.length,sampleNames:sampleNames.slice(0,3),publicVcfFolder,region},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
    // #endregion
    if (!publicVcfFolder || !region || sampleNames.length === 0) {
      // #region agent log
      fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:346',message:'loadRecordsBatch early return',data:{reason:'missing params'},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      return;
    }
    
    // Check which samples need to be loaded
    setLoadedSamples(prev => {
      const toLoad = sampleNames.filter(name => !prev.has(name));
      // #region agent log
      fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:352',message:'loadRecordsBatch check samples',data:{toLoadCount:toLoad.length,prevLoadedCount:prev.size},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      if (toLoad.length === 0) return prev;
      
      // Prevent concurrent loads - set ref synchronously before async work
      if (loadingRecordsRef.current) return prev;
      loadingRecordsRef.current = true;
      setLoadingRecords(true);
      
      // Load the batch asynchronously
      (async () => {
        try {
          const url = `${API_BASE}/api/population/region/${encodeURIComponent(region)}/samples`;
          // #region agent log
          fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:360',message:'loadRecordsBatch API call',data:{url,sample_names:toLoad.join(','),folder_path:publicVcfFolder},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
          // #endregion
          const response = await axios.get(url, {
            params: { 
              folder_path: publicVcfFolder,
              sample_names: toLoad.join(','),
              mode: mode
            }
          });
          // #region agent log
          fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:368',message:'loadRecordsBatch response',data:{status:response.status,records_count:Object.keys(response.data?.records || {}).length},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
          // #endregion
          const records = response.data.records || {};
          setPopulationRecords(prevRecords => ({ ...prevRecords, ...records }));
          setLoadedSamples(prevSet => {
            const newSet = new Set(prevSet);
            toLoad.forEach(name => newSet.add(name));
            return newSet;
          });
        } catch (err: any) {
          // #region agent log
          fetch('http://localhost:7242/ingest/0a9b303e-3b94-470e-983a-030a28b28802',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'CohortAnalysis.tsx:377',message:'loadRecordsBatch error',data:{error:err.message,response_detail:err.response?.data?.detail,status:err.response?.status},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
          // #endregion
          console.error('Error loading records batch:', err);
        } finally {
          loadingRecordsRef.current = false;
          setLoadingRecords(false);
        }
      })();
      
      return prev;
    });
  }, [publicVcfFolder, region, mode]);

  // Load more records progressively
  useEffect(() => {
    if (sampleIds.length === 0 || loadedSamples.size >= sampleIds.length) return;
    // Prevent duplicate loads - check both state and ref (ref is set synchronously)
    if (loadingRecords || loadingRecordsRef.current) return;
    
    // Load next batch if we have unloaded samples
    const unloadedSamples = sampleIds
      .filter(id => !loadedSamples.has(id.sample_name))
      .slice(0, BATCH_SIZE)
      .map(id => id.sample_name);
    
    if (unloadedSamples.length > 0) {
      // Load next batch after a short delay to allow UI to render
      const timer = setTimeout(() => {
        if (!loadingRecordsRef.current) {
          loadRecordsBatch(unloadedSamples);
        }
      }, 500);
      
      return () => clearTimeout(timer);
    }
  }, [sampleIds, loadedSamples, loadingRecords, loadRecordsBatch]);

  // Create a synthetic record object for visualization components
  // Extract motifs from population records
  const syntheticRecord = useMemo(() => {
    if (Object.keys(populationRecords).length === 0) {
      console.log('CohortAnalysis: No population records, syntheticRecord will be null');
      return null;
    }
    
    // Collect all unique motifs from all population records
    const allMotifsSet = new Set<string>();
    Object.values(populationRecords).forEach(rec => {
      if (rec.motifs && Array.isArray(rec.motifs) && rec.motifs.length > 0) {
        rec.motifs.forEach(m => {
          if (m && typeof m === 'string' && m.trim() !== '') {
            allMotifsSet.add(m);
          }
        });
      }
    });
    
    const motifs = Array.from(allMotifsSet);
    
    // Get region info from first record
    const firstRecord = Object.values(populationRecords)[0];
    
    console.log('CohortAnalysis: Creating syntheticRecord', {
      recordCount: Object.keys(populationRecords).length,
      motifCount: motifs.length,
      motifs: motifs.slice(0, 5), // First 5 motifs
      firstRecord: firstRecord ? { chr: firstRecord.chr, pos: firstRecord.pos } : null
    });
    
    return {
      chr: firstRecord?.chr || '',
      pos: firstRecord?.pos || 0,
      stop: firstRecord?.stop || 0,
      motifs: motifs.length > 0 ? motifs : ['Unknown'], // Ensure at least one motif
      motif_ids_h1: [],
      motif_ids_h2: [],
      motif_ids_ref: [],
      ref_CN: null,
      CN_H1: null,
      CN_H2: null,
      spans: [],
      ref_allele: '',
      alt_allele1: '',
      alt_allele2: '',
      gt: '',
      supported_reads_h1: 0,
      supported_reads_h2: 0,
      id: region || 'cohort'
    };
  }, [populationRecords, region]);

  // Get selected sample details
  const selectedSampleData = useMemo(() => {
    if (!selectedSample) return null;
    
    const idx = sequencesData.sequences.findIndex(s => s.name === selectedSample);
    if (idx === -1) return null;
    
    return {
      name: selectedSample,
      sequence: sequencesData.sequences[idx].sequence,
      spans: sequencesData.spanList[idx] || '',
      motifIds: sequencesData.motifIdsList[idx] || [],
      genotype: genotypeComparisonData[selectedSample] || '',
      populationRecord: populationRecords[selectedSample] || 
                        populationRecords[`${selectedSample}_h1`] ||
                        populationRecords[`${selectedSample}_h2`]
    };
  }, [selectedSample, sequencesData, genotypeComparisonData, populationRecords]);

  // Available genotypes
  const availableGenotypes = useMemo(() => {
    const genotypes = Object.values(genotypeComparisonData)
      .filter(gt => gt && gt !== '' && gt !== 'null' && gt !== 'undefined')
      .map(gt => String(gt));
    return Array.from(new Set(genotypes)).sort();
  }, [genotypeComparisonData]);

  // Statistics
  const statistics = useMemo(() => {
    if (sequencesData.sequences.length === 0) {
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
    
    const motifCounts = sequencesData.sequences.map((seq, idx) => {
      const motifIds = sequencesData.motifIdsList[idx] || [];
      return motifIds.filter(id => id && id !== '.' && id !== '').length;
    });
    
    const lengths = sequencesData.sequences.map(seq => seq.sequence.length);
    
    return {
      totalSamples: sequencesData.sortedSamples.length,
      uniqueGenotypes: Object.keys(genotypeGroups).length,
      avgMotifCount: Math.round((motifCounts.reduce((sum, count) => sum + count, 0) / motifCounts.length) * 10) / 10,
      minLength: Math.min(...lengths),
      maxLength: Math.max(...lengths),
      avgLength: Math.round(lengths.reduce((sum, len) => sum + len, 0) / lengths.length)
    };
  }, [sequencesData, genotypeComparisonData]);

  // Genotype frequency
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

  // Handle sample click
  const handleSampleClick = useCallback((sampleName: string) => {
    setSelectedSample(prev => prev === sampleName ? null : sampleName);
  }, []);

  // Extract chromosome and position from region string (format: chr:pos-stop)
  // MUST be called before any early returns (React Hooks rules)
  const regionInfo = useMemo(() => {
    if (!region) return null;
    const match = region.match(/^([\w]+):(\d+)-(\d+)$/);
    if (match) {
      return {
        chr: match[1],
        pos: parseInt(match[2]),
        stop: parseInt(match[3])
      };
    }
    return null;
  }, [region]);
  
  const cohortAnalysisRef = useRef<HTMLDivElement>(null);
  
  const handleExportCSV = useCallback(() => {
    const data = filteredSamples.map(sample => {
      const idx = sequencesData.sequences.findIndex(s => s.name === sample);
      const sequence = idx >= 0 ? sequencesData.sequences[idx] : null;
      const motifIds = idx >= 0 ? sequencesData.motifIdsList[idx] || [] : [];
      const motifCount = motifIds.filter(id => id && id !== '.' && id !== '').length;
      const genotype = genotypeComparisonData[sample] || '';
      
      return {
        Sample: sample,
        Genotype: genotype,
        'Length (bp)': sequence ? sequence.sequence.length : 0,
        'Motif Count': motifCount,
        Sequence: sequence ? sequence.sequence : ''
      };
    });
    
    exportToCSV(data, generateFilename('cohort_analysis', 'csv'));
  }, [filteredSamples, sequencesData, genotypeComparisonData]);
  
  const handleExportFASTA = useCallback(() => {
    const sequences = sequencesData.sequences
      .filter(seq => filteredSamples.includes(seq.name))
      .map(seq => ({
        name: seq.name,
        sequence: seq.sequence
      }));
    exportToFASTA(sequences, generateFilename('cohort_sequences', 'fasta'));
  }, [sequencesData, filteredSamples]);

  if (!publicVcfFolder) {
    return (
      <div className="population-no-data">
        <p>No cohort folder specified. Please load a public VCF folder from the sidebar.</p>
      </div>
    );
  }

  if (error) {
    return <div className="population-error">Error: {error}</div>;
  }

  if (loading && Object.keys(populationRecords).length === 0) {
    return <div className="population-loading">Loading cohort data...</div>;
  }

  if (Object.keys(populationRecords).length === 0) {
    return (
      <div className="population-no-data">
        <p>No cohort data found for region: {region}</p>
        {loading && <p>Still loading...</p>}
      </div>
    );
  }

  // Debug logging
  console.log('CohortAnalysis render:', {
    populationRecordsCount: Object.keys(populationRecords).length,
    sequencesCount: sequencesData.sequences.length,
    syntheticRecordExists: !!syntheticRecord,
    visibleSamplesCount: visibleSamples.length,
    filteredSamplesCount: filteredSamples.length
  });

  // Reuse the CombinedStackHeatmap component from PopulationComparison
  // For now, we'll create a simplified version that makes samples clickable
  
  return (
    <div className="population-comparison" ref={cohortAnalysisRef}>
      <div className="cohort-analysis-header" style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '0.5rem', width: '100%', gap: '1rem' }}>
        <h2 style={{ margin: 0, flex: '1 1 auto', minWidth: 0 }}>👥 Cohort Analysis</h2>
        <div style={{ flexShrink: 0 }}>
          <ExportMenu
            elementRef={cohortAnalysisRef}
            filename="cohort_analysis"
            onExportCSV={handleExportCSV}
            onExportFASTA={handleExportFASTA}
            showImageExport={true}
            showDataExport={true}
          />
        </div>
      </div>
      {regionInfo && (
        <RegionInfoCard
          chr={regionInfo.chr}
          pos={regionInfo.pos}
          stop={regionInfo.stop}
          region={region}
          currentAlleles={[]}
        />
      )}
      {!syntheticRecord && (
        <div className="warning-message" style={{ color: 'orange', margin: '10px 0' }}>
          ⚠️ Warning: Could not create synthetic record for visualizations. Check browser console for details.
        </div>
      )}
      
      {/* Statistics Summary */}
      <div className="population-stats-summary">
        <div className="stat-card-compact">
          <span className="stat-label-compact">Total Samples</span>
          <span className="stat-value-compact">{statistics.totalSamples}</span>
        </div>
        <div className="stat-card-compact">
          <span className="stat-label-compact">Unique Genotypes</span>
          <span className="stat-value-compact">{statistics.uniqueGenotypes}</span>
        </div>
        <div className="stat-card-compact">
          <span className="stat-label-compact">Avg Motif Count</span>
          <span className="stat-value-compact">{statistics.avgMotifCount}</span>
        </div>
        <div className="stat-card-compact">
          <span className="stat-label-compact">Length Range</span>
          <span className="stat-value-compact">{statistics.minLength}-{statistics.maxLength}bp</span>
        </div>
      </div>

      {/* Selected Sample Detail View */}
      {selectedSample && selectedSampleData && (
        <div className="sample-detail-view">
          <div className="sample-detail-header">
            <h3>📋 Sample Details: {selectedSample}</h3>
            <button className="close-btn" onClick={() => setSelectedSample(null)}>✕</button>
          </div>
          <div className="sample-detail-content">
            <div className="sample-detail-info">
              <div><strong>Genotype:</strong> {selectedSampleData.genotype || 'N/A'}</div>
              <div><strong>Sequence Length:</strong> {selectedSampleData.sequence.length} bp</div>
              <div><strong>Motif Count:</strong> {selectedSampleData.motifIds.filter(id => id && id !== '.' && id !== '').length}</div>
              {selectedSampleData.populationRecord && (
                <>
                  <div><strong>Reference CN:</strong> {selectedSampleData.populationRecord.ref_CN}</div>
                  <div><strong>Haplotype CN:</strong> {selectedSampleData.populationRecord.CN_H}</div>
                  <div><strong>Region ID:</strong> {selectedSampleData.populationRecord.id}</div>
                </>
              )}
            </div>
            <div className="sample-sequence-view">
              <h4>Sequence:</h4>
              <div className="sequence-display">
                {selectedSampleData.sequence}
              </div>
            </div>
          </div>
        </div>
      )}

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

      {/* Genotype Frequency Distribution */}
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
                    />
                  </div>
                </div>
              ))}
            </div>
          </div>
        )}
      </div>

      {/* Search and Filter */}
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
      {syntheticRecord && (
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
                  Showing {visibleSamples.length} of {filteredSamples.length} samples
                  <div className="virtual-scroll-controls">
                    <button
                      className="virtual-scroll-btn"
                      onClick={() => {
                        const populationCount = filteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length;
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
                        const populationCount = filteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length;
                        setVisibleRange(prev => ({
                          start: prev.start + 20,
                          end: Math.min(populationCount, prev.end + 20)
                        }));
                      }}
                      disabled={visibleRange.end >= filteredSamples.filter(s => s !== 'Ref' && !s.startsWith('Allel')).length}
                    >
                      ↓ Load More Below
                    </button>
                  </div>
                </div>
              )}
              {syntheticRecord && sequencesData.sequences.length > 0 && visibleSamples.length > 0 ? (
                (() => {
                  // Filter sequences and keep corresponding spanList and motifIdsList aligned
                  const filteredSequences: Array<{ name: string; sequence: string }> = [];
                  const filteredSpanList: string[] = [];
                  const filteredMotifIdsList: string[][] = [];
                  
                  console.log('CohortAnalysis: Filtering sequences for plot', {
                    totalSequences: sequencesData.sequences.length,
                    visibleSamplesCount: visibleSamples.length,
                    sequenceNames: sequencesData.sequences.map(s => s.name).slice(0, 5),
                    visibleSampleNames: visibleSamples.slice(0, 5)
                  });
                  
                  // Helper function to extract base name (remove _h1, _h2, etc.)
                  const getBaseName = (name: string): string => {
                    // Remove _h1, _h2, _hap1, _hap2, _haplotype_1, _haplotype_2, etc.
                    return name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
                  };
                  
                  // Filter sequences based on filteredSamples (which is already sorted correctly)
                  sequencesData.sequences.forEach((seq, idx) => {
                    const baseName = getBaseName(seq.name);
                    // Include if exact match or base name matches (also always include Ref)
                    if (seq.name === 'Ref' || filteredSamples.includes(seq.name) || filteredSamples.includes(baseName)) {
                      filteredSequences.push(seq);
                      filteredSpanList.push(sequencesData.spanList[idx] || '');
                      filteredMotifIdsList.push(sequencesData.motifIdsList[idx] || []);
                    }
                  });
                  
                  // Sort filteredSequences to match the order of filteredSamples (which respects user's sort selection)
                  const sortedFilteredSequences = [...filteredSequences].sort((a, b) => {
                    // Get base names for matching
                    const baseNameA = getBaseName(a.name);
                    const baseNameB = getBaseName(b.name);
                    
                    // Find index in filteredSamples (try exact match first, then base name)
                    let indexA = filteredSamples.indexOf(a.name);
                    if (indexA === -1) {
                      indexA = filteredSamples.indexOf(baseNameA);
                    }
                    
                    let indexB = filteredSamples.indexOf(b.name);
                    if (indexB === -1) {
                      indexB = filteredSamples.indexOf(baseNameB);
                    }
                    
                    // If still not found, put at end
                    if (indexA === -1 && indexB === -1) return 0;
                    if (indexA === -1) return 1;
                    if (indexB === -1) return -1;
                    
                    // If same base name position, sort by haplotype (_h1 before _h2) to maintain consistent order
                    if (indexA === indexB && baseNameA === baseNameB) {
                      const hasH1A = a.name.includes('_h1') || a.name.includes('_hap1');
                      const hasH1B = b.name.includes('_h1') || b.name.includes('_hap1');
                      if (hasH1A && !hasH1B) return -1;
                      if (!hasH1A && hasH1B) return 1;
                      // If both have same haplotype suffix or neither, maintain original order
                      return a.name.localeCompare(b.name);
                    }
                    
                    return indexA - indexB;
                  });
                  
                  // Reorder spanList and motifIdsList to match sorted sequences
                  const sortedFilteredSpanList: string[] = [];
                  const sortedFilteredMotifIdsList: string[][] = [];
                  sortedFilteredSequences.forEach(seq => {
                    const originalIdx = sequencesData.sequences.findIndex(s => s.name === seq.name);
                    if (originalIdx >= 0) {
                      sortedFilteredSpanList.push(sequencesData.spanList[originalIdx] || '');
                      sortedFilteredMotifIdsList.push(sequencesData.motifIdsList[originalIdx] || []);
                    }
                  });
                  
                  // Use the sorted sequence names to maintain proper order
                  const matchingSequenceNames = sortedFilteredSequences.map(seq => seq.name);
                  
                  console.log('CohortAnalysis: After filtering and sorting', {
                    filteredSequencesCount: sortedFilteredSequences.length,
                    matchingSequenceNamesCount: matchingSequenceNames.length,
                    matchingSequenceNames: matchingSequenceNames.slice(0, 5),
                    sortBy,
                    sortDirection
                  });
                  
                  if (sortedFilteredSequences.length === 0) {
                    return (
                      <div className="info-message" style={{ padding: '20px', textAlign: 'center' }}>
                        <p>⚠️ No sequences match visible samples after filtering</p>
                        <p style={{ fontSize: '0.9em', marginTop: '10px', color: '#888' }}>
                          Debug: {sequencesData.sequences.length} total sequences, {visibleSamples.length} visible samples
                        </p>
                        <p style={{ fontSize: '0.8em', marginTop: '5px', color: '#666' }}>
                          Sequence names: {sequencesData.sequences.map(s => s.name).slice(0, 10).join(', ')}
                        </p>
                        <p style={{ fontSize: '0.8em', marginTop: '5px', color: '#666' }}>
                          Visible samples: {visibleSamples.slice(0, 10).join(', ')}
                        </p>
                      </div>
                    );
                  }
                  
                  return (
                    <CombinedStackHeatmap 
                      record={syntheticRecord}
                      sequences={sortedFilteredSequences}
                      spanList={sortedFilteredSpanList}
                      motifIdsList={sortedFilteredMotifIdsList}
                      sortedSamples={matchingSequenceNames}
                      pathogenicThreshold={pathogenicInfo?.pathogenic_threshold}
                    />
                  );
                })()
              ) : (
                <div className="info-message" style={{ padding: '20px', textAlign: 'center' }}>
                  {!syntheticRecord && <p>⚠️ No synthetic record available. Check console for details.</p>}
                  {syntheticRecord && sequencesData.sequences.length === 0 && <p>⚠️ No sequences extracted from population records. Check data structure.</p>}
                  {syntheticRecord && sequencesData.sequences.length > 0 && visibleSamples.length === 0 && <p>⚠️ No visible samples after filtering.</p>}
                  <p style={{ fontSize: '0.9em', marginTop: '10px', color: '#888' }}>
                    Debug: {sequencesData.sequences.length} sequences, {visibleSamples.length} visible samples
                  </p>
                </div>
              )}
            </>
          )}
        </div>
      )}

      {/* Bar Plot */}
      {syntheticRecord && (
        <div className="population-section">
          <div 
            className="section-header-collapsible"
            onClick={() => setExpandedSections({...expandedSections, barPlot: !expandedSections.barPlot})}
          >
            <h3>Motif Count per Sample</h3>
            <span className="collapse-icon">{expandedSections.barPlot ? '▼' : '▶'}</span>
          </div>
          {expandedSections.barPlot && (
            sequencesData.sequences.length > 0 ? (
              <BarPlot 
                sequences={sequencesData.sequences} 
                motifIdsList={sequencesData.motifIdsList}
                sortedSamples={filteredSamples.filter(s => 
                  sequencesData.sequences.some(seq => seq.name === s)
                )}
                pathogenicThreshold={pathogenicInfo?.pathogenic_threshold}
              />
            ) : (
              <div className="info-message">No sequences available for bar plot</div>
            )
          )}
        </div>
      )}

      {/* Cluster Plot */}
      {syntheticRecord && (
        <div className="population-section">
          <div 
            className="section-header-collapsible"
            onClick={() => setExpandedSections({...expandedSections, cluster: !expandedSections.cluster})}
          >
            <h3>Sample Clustering</h3>
            <span className="collapse-icon">{expandedSections.cluster ? '▼' : '▶'}</span>
          </div>
          {expandedSections.cluster && (
            sequencesData.sequences.length > 0 ? (
              <ClusterPlot sequences={sequencesData.sequences} motifIdsList={sequencesData.motifIdsList} />
            ) : (
              <div className="info-message">No sequences available for cluster plot</div>
            )
          )}
        </div>
      )}

      {/* Sample Table with Sortable Columns and Quick Filters */}
      <div className="population-section">
        <div 
          className="section-header-collapsible"
          onClick={() => setExpandedSections({...expandedSections, samples: !expandedSections.samples})}
        >
          <h3>Cohort Samples ({filteredSamples.length})</h3>
          <span className="collapse-icon">{expandedSections.samples ? '▼' : '▶'}</span>
        </div>
        {expandedSections.samples && (
          <>
            {/* Quick Filters for table */}
            <div className="cohort-table-filters">
              <div className="cohort-filter-group">
                <span className="cohort-filter-label">Status:</span>
                <div className="cohort-filter-buttons">
                  {(['all', 'loaded', 'not-loaded'] as const).map(f => (
                    <button
                      key={f}
                      className={`cohort-filter-btn ${statusFilter === f ? 'active' : ''}`}
                      onClick={() => setStatusFilter(f)}
                    >
                      {f === 'all' ? 'All' : f === 'loaded' ? 'Loaded' : 'Not Loaded'}
                    </button>
                  ))}
                </div>
              </div>
              <div className="cohort-filter-group">
                <span className="cohort-filter-label">Genotype:</span>
                <div className="cohort-filter-buttons">
                  {availableGenotypes.map(gt => (
                    <button
                      key={gt}
                      className={`cohort-filter-btn small ${selectedGenotypes.includes(gt) ? 'active' : ''}`}
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
                      className="cohort-filter-btn small clear"
                      onClick={() => setSelectedGenotypes([])}
                    >
                      Clear
                    </button>
                  )}
                </div>
              </div>
              <div className="cohort-filter-group cohort-search-inline">
                <input
                  type="text"
                  className="cohort-table-search"
                  placeholder="Search samples..."
                  value={searchQuery}
                  onChange={(e) => setSearchQuery(e.target.value)}
                />
              </div>
            </div>
            {/* Sortable Table */}
            <div className="cohort-table-wrapper">
              <table className="cohort-samples-table">
                <thead>
                  <tr>
                    <th
                      className="cohort-table-th sortable"
                      onClick={() => {
                        if (sortBy === 'name') {
                          setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                        } else {
                          setSortBy('name');
                          setSortDirection('asc');
                        }
                      }}
                      title="Click to sort"
                    >
                      Sample {sortBy === 'name' && (sortDirection === 'asc' ? '↑' : '↓')}
                    </th>
                    <th
                      className="cohort-table-th sortable"
                      onClick={() => {
                        if (sortBy === 'genotype') {
                          setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                        } else {
                          setSortBy('genotype');
                          setSortDirection('asc');
                        }
                      }}
                      title="Click to sort"
                    >
                      Genotype {sortBy === 'genotype' && (sortDirection === 'asc' ? '↑' : '↓')}
                    </th>
                    <th
                      className="cohort-table-th sortable"
                      onClick={() => {
                        if (sortBy === 'length') {
                          setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                        } else {
                          setSortBy('length');
                          setSortDirection('desc');
                        }
                      }}
                      title="Click to sort"
                    >
                      Length (bp) {sortBy === 'length' && (sortDirection === 'asc' ? '↑' : '↓')}
                    </th>
                    <th
                      className="cohort-table-th sortable"
                      onClick={() => {
                        if (sortBy === 'motifCount') {
                          setSortDirection(sortDirection === 'asc' ? 'desc' : 'asc');
                        } else {
                          setSortBy('motifCount');
                          setSortDirection('desc');
                        }
                      }}
                      title="Click to sort"
                    >
                      Motifs {sortBy === 'motifCount' && (sortDirection === 'asc' ? '↑' : '↓')}
                    </th>
                    <th className="cohort-table-th">Status</th>
                  </tr>
                </thead>
                <tbody>
                  {filteredSamples.map(sample => {
                    const isLoaded = loadedSamples.has(sample) || sample === 'Ref';
                    // Find sequence data: exact match, or base-name match for diploid (e.g. NA12878 -> NA12878_h1, NA12878_h2)
                    const getBaseName = (name: string) => name.replace(/_(h[12]|hap[12]|haplotype_[12])$/i, '');
                    const sampleBase = getBaseName(sample);
                    const matchingIndices = sequencesData.sequences
                      .map((s, i) => ({ seq: s, i }))
                      .filter(({ seq }) => seq.name === sample || getBaseName(seq.name) === sampleBase);
                    let totalLength = 0;
                    let motifCount = 0;
                    if (matchingIndices.length > 0) {
                      totalLength = Math.max(...matchingIndices.map(({ seq }) => seq.sequence?.length || 0));
                      motifCount = matchingIndices.reduce((sum, { i }) => {
                        const ids = sequencesData.motifIdsList[i] || [];
                        return sum + ids.filter(id => id && id !== '.' && id !== '').length;
                      }, 0);
                    }
                    const genotype = genotypeComparisonData[sample] || '';
                    const isSelected = selectedSample === sample;
                    const isLoading = !isLoaded && loadingRecords;
                    const isRef = sample === 'Ref' || sample.startsWith('Allel');
                    
                    return (
                      <tr
                        key={sample}
                        className={`cohort-table-row ${isSelected ? 'selected' : ''} ${isRef ? 'reference' : ''} ${!isLoaded ? 'not-loaded' : ''}`}
                        onClick={() => {
                          if (!isLoaded && !isLoading) {
                            loadRecordsBatch([sample]);
                          }
                          handleSampleClick(sample);
                        }}
                        title={isLoaded 
                          ? `Click to view details: ${sample}`
                          : `Click to load: ${sample}`
                        }
                      >
                        <td className="cohort-table-td sample-name">
                          {sample}
                          {!isLoaded && <span className="loading-badge"> ⏳</span>}
                        </td>
                        <td className="cohort-table-td">{genotype || 'N/A'}</td>
                        <td className="cohort-table-td numeric">
                          {totalLength > 0 ? totalLength : '—'}
                        </td>
                        <td className="cohort-table-td numeric">{motifCount}</td>
                        <td className="cohort-table-td">
                          <span className={`status-badge ${isLoaded ? 'loaded' : 'pending'}`}>
                            {isLoaded ? 'Loaded' : 'Click to load'}
                          </span>
                        </td>
                      </tr>
                    );
                  })}
                </tbody>
              </table>
            </div>
            {loadingRecords && (
              <div className="cohort-loading-more">
                Loading more samples... ({loadedSamples.size} / {sampleIds.length} loaded)
              </div>
            )}
          </>
        )}
      </div>
    </div>
  );
};

// Helper functions
function computeDiploidGenotype(gt_h1: string, gt_h2: string, ids_h1: string[], ids_h2: string[]): string {
  const g1 = parseInt(gt_h1) || 0;
  const g2 = parseInt(gt_h2) || 0;
  
  if (g1 === 0 && g2 === 0) return '0/0';
  if ((g1 === 0 && g2 === 1) || (g1 === 1 && g2 === 0)) return '0/1';
  if (g1 === 1 && g2 === 1) {
    const ids1Str = ids_h1.join('_');
    const ids2Str = ids_h2.join('_');
    return ids1Str === ids2Str ? '1/1' : '1/2';
  }
  return `${g1}/${g2}`;
}

function parseMotifRange(spanStr: string): Array<[number, number]> {
  if (!spanStr) return [];
  
  // Check cache first
  if (parsedRangeCache.has(spanStr)) {
    return parsedRangeCache.get(spanStr)!;
  }
  
  const pattern = /\((\d+)-(\d+)\)/g;
  const ranges: Array<[number, number]> = [];
  let match: RegExpExecArray | null;
  while ((match = pattern.exec(spanStr)) !== null) {
    ranges.push([parseInt(match[1]) - 1, parseInt(match[2]) - 1]);
  }
  
  // Cache the result (limit cache size)
  if (parsedRangeCache.size > 1000) {
    const firstKey = parsedRangeCache.keys().next().value;
    if (firstKey) parsedRangeCache.delete(firstKey);
  }
  parsedRangeCache.set(spanStr, ranges);
  
  return ranges;
}

// Interface for Record used by visualization components
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

// Combined Stack Plot and Heatmap Component
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
  const plotCacheKey = useMemo(() => {
    if (!sequences || sequences.length === 0) {
      return null; // Don't create cache key if no sequences
    }
    const seqKey = sequences.map(s => s.name).join(',');
    const spanKey = spanList.join('|');
    const motifKey = motifIdsList.map(m => m.join('_')).join('|');
    return `${record.id}:${seqKey}:${spanKey}:${motifKey}`;
  }, [record.id, sequences, spanList, motifIdsList]);
  
  const plotData = useMemo(() => {
    // Don't compute plot data if we don't have sequences
    if (!sequences || sequences.length === 0) {
      console.log('CohortAnalysis: Skipping plot data computation - no sequences', {
        sequencesCount: sequences?.length || 0
      });
      return null;
    }
    
    if (plotCacheKey) {
      const cached = combinedPlotCache.get(plotCacheKey);
      if (cached) {
        console.log('CohortAnalysis: Using cached plot data', { 
          motifDataCount: cached.motifData.length, 
          heatmapDataCount: cached.heatmapData.length 
        });
        return cached;
      }
    }
    
    console.log('CohortAnalysis: Computing plot data', {
      sequencesCount: sequences.length,
      spanListCount: spanList.length,
      motifIdsListCount: motifIdsList.length,
      recordMotifs: record.motifs?.length || 0,
      sortedSamplesCount: sortedSamples.length
    });
    
    const stackData: Array<{
      sample: string;
      start: number;
      end: number;
      motif: string;
      motifIndex: number;
      length: number;
    }> = [];
    
    for (let seqIdx = 0; seqIdx < sequences.length; seqIdx++) {
      const seq = sequences[seqIdx];
      const span = spanList[seqIdx] || '';
      const motifIds = motifIdsList[seqIdx] || [];
      const ranges = parseMotifRange(span);
      
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
        const motifName = (record.motifs && record.motifs.length > 0 && record.motifs[motifId]) 
          ? record.motifs[motifId] 
          : 'Unknown';
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
    
    const heatData: Array<{ sample: string; motif: string; count: number }> = [];
    
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
  
    const motifTotalsMap = new Map<string, number>();
    for (const d of heatData) {
      motifTotalsMap.set(d.motif, (motifTotalsMap.get(d.motif) || 0) + d.count);
    }
    
    const motifTotals = Array.from(motifTotalsMap.entries())
      .map(([motif, total]) => ({ motif, total }))
      .sort((a, b) => b.total - a.total);
    const sortedMots = motifTotals.map(m => m.motif);
    
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
    
    // Calculate maxLen from actual sequence lengths, not from stackData (which is based on SP ranges)
    let maxLen = 1;
    for (const seq of sequences) {
      if (seq.sequence && seq.sequence.length > maxLen) {
        maxLen = seq.sequence.length;
      }
    }
    
    let maxCnt = 1;
    for (const d of heatData) {
      if (d.count > maxCnt) maxCnt = d.count;
    }
    
    const sampleDataMap: { [sample: string]: Array<{
      start: number;
      end: number;
      motif: string;
      motifIndex: number;
      length: number;
    }> } = {};
    
    for (const sample of sortedSamples) {
      sampleDataMap[sample] = [];
    }
    
    for (const item of stackData) {
      if (sampleDataMap[item.sample]) {
        sampleDataMap[item.sample].push(item);
      }
    }
    
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
    
    if (plotCacheKey) {
      combinedPlotCache.set(plotCacheKey, result);
      
      if (combinedPlotCache.size > 50) {
        const firstKey = combinedPlotCache.keys().next().value;
        if (firstKey) {
          combinedPlotCache.delete(firstKey);
        }
      }
    }
    
    return result;
  }, [sequences, spanList, motifIdsList, record.motifs, plotCacheKey, sortedSamples]);
  
  const { motifData = [], heatmapData = [], motifColorMap = {}, maxLength = 1, maxCount = 1, sortedMotifs = [], sampleDataMap = {} } = plotData || {};
  
  const heatmapDataMap = useMemo(() => {
    const map: { [key: string]: number } = {};
    heatmapData.forEach(d => {
      map[`${d.sample}_${d.motif}`] = d.count;
    });
    return map;
  }, [heatmapData]);
  
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
  
  const uniqueMotifsForLegend = useMemo(() => {
    const uniqueMotifs = new Set<string>();
    for (const d of motifData) {
      if (d.motif !== 'Interruption') {
        uniqueMotifs.add(d.motif);
      }
    }
    return Array.from(uniqueMotifs);
  }, [motifData]);
  
  // Early return if we don't have the required data - check sequences first
  if (!sequences || sequences.length === 0) {
    console.log('CohortAnalysis: CombinedStackHeatmap early return - no sequences', {
      sequencesCount: sequences?.length || 0
    });
    return (
      <div className="info-message">
        <p>No sequences available to display</p>
        <p style={{ fontSize: '0.9em', marginTop: '10px', color: '#888' }}>
          Debug: sequences={sequences?.length || 0}
        </p>
      </div>
    );
  }
  
  if (!sortedSamples || sortedSamples.length === 0) {
    console.log('CohortAnalysis: CombinedStackHeatmap early return - no sorted samples', {
      sortedSamplesCount: sortedSamples?.length || 0,
      sequencesCount: sequences.length
    });
    return (
      <div className="info-message">
        <p>No samples available to display</p>
        <p style={{ fontSize: '0.9em', marginTop: '10px', color: '#888' }}>
          Debug: sortedSamples={sortedSamples?.length || 0}, sequences={sequences.length}
        </p>
      </div>
    );
  }
  
  if (!plotData) {
    console.log('CohortAnalysis: CombinedStackHeatmap early return - no plot data', {
      sequencesCount: sequences.length,
      sortedSamplesCount: sortedSamples.length
    });
    return (
      <div className="info-message">
        <p>Computing plot data...</p>
        <p style={{ fontSize: '0.9em', marginTop: '10px', color: '#888' }}>
          Debug: sequences={sequences.length}, sortedSamples={sortedSamples.length}
        </p>
      </div>
    );
  }
  
  return (
    <div className="combined-stack-heatmap">
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
          const sampleStackData = sampleDataMap[sample] || [];
          // Use actual sequence length instead of calculating from stack data (which only includes motif spans)
          const sequenceIdx = sequences.findIndex(s => s.name === sample);
          const actualSequenceLength = sequenceIdx >= 0 && sequences[sequenceIdx] ? sequences[sequenceIdx].sequence.length : 0;
          // Calculate displayed length from stack data for visualization scaling
          let displayedLength = 0;
          for (const d of sampleStackData) {
            if (d.end > displayedLength) displayedLength = d.end;
          }
          // Use actual sequence length if available, otherwise fall back to displayed length
          const totalLength = actualSequenceLength > 0 ? actualSequenceLength : (displayedLength > 0 ? displayedLength : 1);
          
          // Check if this sample exceeds pathogenic threshold (using motif count, not sequence length)
          const exceedsThreshold = pathogenicThreshold !== undefined && pathogenicThreshold !== null && pathogenicThreshold > 0 && (() => {
            const thresholdNum = typeof pathogenicThreshold === 'number' ? pathogenicThreshold : parseFloat(String(pathogenicThreshold));
            if (isNaN(thresholdNum) || thresholdNum <= 0) return false;
            
            // Get motif count for this sample
            const motifIds = sequenceIdx >= 0 && motifIdsList[sequenceIdx] ? motifIdsList[sequenceIdx] : [];
            const motifCount = motifIds.filter(id => id && id !== '.' && id !== '').length;
            
            return motifCount >= thresholdNum;
          })();
          
          return (
            <div 
              key={sample} 
              className={`combined-row ${sample === 'Ref' || sample.startsWith('Allel') ? 'current-sample-row' : ''} ${exceedsThreshold ? 'pathogenic-row' : ''}`}
              style={exceedsThreshold ? {
                border: '3px solid #DC2626',
                borderRadius: '4px',
                boxShadow: '0 0 8px rgba(220, 38, 38, 0.4)'
              } : {}}
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
                  
                  return (
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

// Bar Plot Component
interface BarPlotProps {
  sequences: Array<{ name: string; sequence: string }>;
  motifIdsList: string[][];
  sortedSamples?: string[];
  pathogenicThreshold?: number;
}

const BarPlot: React.FC<BarPlotProps> = memo(({ sequences, motifIdsList, sortedSamples, pathogenicThreshold }) => {
  // Create a map for quick lookup
  const sequenceMap = useMemo(() => {
    const map = new Map<string, { sequence: string; index: number }>();
    sequences.forEach((seq, idx) => {
      map.set(seq.name, { sequence: seq.sequence, index: idx });
    });
    return map;
  }, [sequences]);

  const motifCounts = useMemo(() => {
    // Use sortedSamples if provided, otherwise use sequences order (maintains consistent sorting across plots)
    const sampleOrder = sortedSamples && sortedSamples.length > 0 
      ? sortedSamples 
      : sequences.map(s => s.name);
    
    return sampleOrder.map((sampleName) => {
      const seqData = sequenceMap.get(sampleName);
      if (!seqData) return { sample: sampleName, count: 0 };
      
      const motifIds = motifIdsList[seqData.index] || [];
      const validMotifs = motifIds.filter(id => id && id !== '.' && id !== '');
      return {
        sample: sampleName,
        count: validMotifs.length
      };
    });
  }, [sequences, motifIdsList, sortedSamples, sequenceMap]);
  
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
  
  const stats = useMemo(() => {
    const total = motifCounts.reduce((sum, item) => sum + item.count, 0);
    const avg = total / motifCounts.length;
    const sorted = [...motifCounts].sort((a, b) => a.count - b.count);
    const median = sorted.length > 0 
      ? sorted[Math.floor(sorted.length / 2)].count 
      : 0;
    return { total, avg: Math.round(avg * 10) / 10, median };
  }, [motifCounts]);
  
  const yAxisTicks = useMemo(() => {
    const numTicks = 6;
    const ticks: number[] = [];
    for (let i = 0; i <= numTicks; i++) {
      const value = Math.round((effectiveMaxCount / numTicks) * i);
      ticks.push(value);
    }
    return ticks.reverse();
  }, [effectiveMaxCount]);
  
  const CHART_HEIGHT = 350;
  const BAR_WIDTH = 60;
  const BAR_GAP = 20;
  const Y_AXIS_WIDTH = 60;
  
  return (
    <div className="bar-plot-v2">
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
      
      <div className="bar-plot-chart-new">
        <div className="chart-unified-container" style={{ height: `${CHART_HEIGHT}px` }}>
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
          
          <div className="chart-area-new" style={{ flex: 1, height: `${CHART_HEIGHT}px`, position: 'relative' }}>
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

// K-Means Clustering with Schuilliot Method (Silhouette Score)
// Helper functions for k-means clustering
function standardScaler(data: number[][]): { scaled: number[][]; mean: number[]; std: number[] } {
  if (data.length === 0) return { scaled: [], mean: [], std: [] };
  
  const n = data[0].length;
  const mean: number[] = [];
  const std: number[] = [];
  
  // Calculate mean for each dimension
  for (let i = 0; i < n; i++) {
    const values = data.map(row => row[i]);
    mean[i] = values.reduce((a, b) => a + b, 0) / values.length;
  }
  
  // Calculate standard deviation for each dimension
  for (let i = 0; i < n; i++) {
    const values = data.map(row => row[i]);
    const variance = values.reduce((sum, val) => sum + Math.pow(val - mean[i], 2), 0) / values.length;
    std[i] = Math.sqrt(variance) || 1; // Avoid division by zero
  }
  
  // Scale the data
  const scaled = data.map(row => 
    row.map((val, i) => (val - mean[i]) / std[i])
  );
  
  return { scaled, mean, std };
}

function euclideanDistance(a: number[], b: number[]): number {
  return Math.sqrt(a.reduce((sum, val, i) => sum + Math.pow(val - b[i], 2), 0));
}

function kMeans(data: number[][], k: number, maxIterations: number = 100, randomSeed: number = 42): number[] {
  if (data.length === 0 || k <= 0 || k > data.length) return [];
  
  // Simple seeded random for reproducibility
  let seed = randomSeed;
  const random = () => {
    seed = (seed * 9301 + 49297) % 233280;
    return seed / 233280;
  };
  
  // Initialize centroids randomly
  const centroids: number[][] = [];
  const usedIndices = new Set<number>();
  for (let i = 0; i < k; i++) {
    let idx;
    do {
      idx = Math.floor(random() * data.length);
    } while (usedIndices.has(idx));
    usedIndices.add(idx);
    centroids.push([...data[idx]]);
  }
  
  let labels: number[] = [];
  let changed = true;
  let iterations = 0;
  
  while (changed && iterations < maxIterations) {
    // Assign points to nearest centroid
    const newLabels = data.map(point => {
      let minDist = Infinity;
      let nearestCluster = 0;
      centroids.forEach((centroid, idx) => {
        const dist = euclideanDistance(point, centroid);
        if (dist < minDist) {
          minDist = dist;
          nearestCluster = idx;
        }
      });
      return nearestCluster;
    });
    
    // Check if labels changed
    changed = !labels.every((label, i) => label === newLabels[i]);
    labels = newLabels;
    
    // Update centroids
    for (let i = 0; i < k; i++) {
      const clusterPoints = data.filter((_, idx) => labels[idx] === i);
      if (clusterPoints.length > 0) {
        const n = data[0].length;
        centroids[i] = [];
        for (let j = 0; j < n; j++) {
          const sum = clusterPoints.reduce((acc, point) => acc + point[j], 0);
          centroids[i][j] = sum / clusterPoints.length;
        }
      }
    }
    
    iterations++;
  }
  
  return labels;
}

function silhouetteScore(data: number[][], labels: number[]): number {
  if (data.length === 0 || labels.length === 0) return -1;
  
  const n = data.length;
  const uniqueLabels = Array.from(new Set(labels));
  if (uniqueLabels.length < 2) return -1; // Need at least 2 clusters
  
  let totalScore = 0;
  
  for (let i = 0; i < n; i++) {
    const point = data[i];
    const label = labels[i];
    
    // Calculate average distance to points in same cluster (a_i)
    const sameClusterPoints = data.filter((_, idx) => labels[idx] === label && idx !== i);
    const a_i = sameClusterPoints.length > 0
      ? sameClusterPoints.reduce((sum, p) => sum + euclideanDistance(point, p), 0) / sameClusterPoints.length
      : 0;
    
    // Calculate minimum average distance to points in other clusters (b_i)
    let min_b_i = Infinity;
    for (const otherLabel of uniqueLabels) {
      if (otherLabel === label) continue;
      
      const otherClusterPoints = data.filter((_, idx) => labels[idx] === otherLabel);
      if (otherClusterPoints.length > 0) {
        const b_i = otherClusterPoints.reduce((sum, p) => sum + euclideanDistance(point, p), 0) / otherClusterPoints.length;
        min_b_i = Math.min(min_b_i, b_i);
      }
    }
    
    if (min_b_i === Infinity) min_b_i = a_i; // Fallback
    
    // Silhouette score for this point
    const s_i = min_b_i > a_i ? (min_b_i - a_i) / Math.max(a_i, min_b_i) : 0;
    totalScore += s_i;
  }
  
  return totalScore / n;
}

function calculateInertia(data: number[][], labels: number[], centroids: number[][]): number {
  let inertia = 0;
  data.forEach((point, idx) => {
    const centroid = centroids[labels[idx]];
    if (centroid) {
      inertia += Math.pow(euclideanDistance(point, centroid), 2);
    }
  });
  return inertia;
}

// Cluster Plot Component with K-Means (Schuilliot Method)
interface ClusterPlotProps {
  sequences: Array<{ name: string; sequence: string }>;
  motifIdsList: string[][];
}

const ClusterPlot: React.FC<ClusterPlotProps> = memo(({ sequences, motifIdsList }) => {
  const [clusterBy, setClusterBy] = useState<'Copy Number' | 'Length' | 'Both'>('Both');
  const [numClusters, setNumClusters] = useState<number | null>(null); // null means use suggested
  const [userChangedK, setUserChangedK] = useState<boolean>(false); // Track if user manually changed K
  
  // Calculate sample statistics
  const sampleStats = useMemo(() => {
    return sequences
      .map((seq, idx) => {
        const motifIds = motifIdsList[idx] || [];
        const validMotifs = motifIds.filter(id => id && id !== '.' && id !== '');
        const copyNumber = validMotifs.length;
        const length = seq.sequence.length;
        
        return {
          sample: seq.name,
          copyNumber,
          length
        };
      })
      .filter(stat => stat.sample !== 'Interruption');
  }, [sequences, motifIdsList]);
  
  // Prepare data for clustering
  const { clusterData, scaledData, suggestedK, metrics } = useMemo(() => {
    if (sampleStats.length < 2) {
      return { clusterData: [], scaledData: [], suggestedK: 2, metrics: [] };
    }
    
    let X: number[][];
    let X_use: number[][];
    
    if (clusterBy === 'Both') {
      X = sampleStats.map(s => [s.copyNumber, s.length]);
      const scaled = standardScaler(X);
      X_use = scaled.scaled;
    } else if (clusterBy === 'Copy Number') {
      X = sampleStats.map(s => [s.copyNumber]);
      X_use = X;
    } else {
      X = sampleStats.map(s => [s.length]);
      X_use = X;
    }
    
    // Test different k values for silhouette score
    const minK = 2;
    const maxK = Math.min(8, sampleStats.length);
    const metrics: Array<{ k: number; actualK: number; silhouette: number; inertia: number }> = [];
    
    for (let k = minK; k <= maxK; k++) {
      if (k >= sampleStats.length) break;
      
      const labels = kMeans(X_use, k, 100, 42);
      const actualK = Array.from(new Set(labels)).length; // Count actual clusters created
      const silhouette = silhouetteScore(X_use, labels);
      
      // Calculate inertia
      const centroids: number[][] = [];
      for (let i = 0; i < k; i++) {
        const clusterPoints = X_use.filter((_, idx) => labels[idx] === i);
        if (clusterPoints.length > 0) {
          const n = X_use[0].length;
          centroids[i] = [];
          for (let j = 0; j < n; j++) {
            const sum = clusterPoints.reduce((acc, point) => acc + point[j], 0);
            centroids[i][j] = sum / clusterPoints.length;
          }
        }
      }
      const inertia = calculateInertia(X_use, labels, centroids);
      
      metrics.push({ k, actualK, silhouette, inertia });
    }
    
    // Suggest k with max silhouette score, but prefer metrics where actualK == k (no empty clusters)
    const validMetrics = metrics.filter(m => !isNaN(m.silhouette) && m.silhouette !== -1);
    let suggestedK = 2;
    if (validMetrics.length > 0) {
      // Prefer metrics where actualK == k (no empty clusters)
      const perfectMetrics = validMetrics.filter(m => m.actualK === m.k);
      const metricsToUse = perfectMetrics.length > 0 ? perfectMetrics : validMetrics;
      
      const maxSilhouette = Math.max(...metricsToUse.map(m => m.silhouette));
      const bestMetric = metricsToUse.find(m => m.silhouette === maxSilhouette);
      if (bestMetric) {
        suggestedK = bestMetric.k;
      }
    } else {
      // Fallback to elbow method
      if (metrics.length > 1) {
        const diffs = [];
        for (let i = 1; i < metrics.length; i++) {
          diffs.push(metrics[i - 1].inertia - metrics[i].inertia);
        }
        if (diffs.length > 0) {
          const maxDiff = Math.max(...diffs);
          const elbowIdx = diffs.indexOf(maxDiff);
          suggestedK = metrics[elbowIdx + 1]?.k || 2;
        }
      }
    }
    
    return {
      clusterData: X,
      scaledData: X_use,
      suggestedK,
      metrics
    };
  }, [sampleStats, clusterBy]);
  
  // Determine the actual K to use: user's choice if they changed it, otherwise suggestedK
  const actualK = useMemo(() => {
    if (numClusters !== null && userChangedK) {
      return numClusters;
    }
    return suggestedK;
  }, [numClusters, userChangedK, suggestedK]);
  
  // Update numClusters when suggestedK changes (only when clusterBy changes, reset user choice)
  const prevClusterBy = useRef(clusterBy);
  useEffect(() => {
    // Reset to suggested K when clusterBy changes (user is changing the clustering method)
    if (prevClusterBy.current !== clusterBy) {
      setNumClusters(null); // Will use suggestedK
      setUserChangedK(false);
      prevClusterBy.current = clusterBy;
    }
  }, [clusterBy]);
  
  // Perform clustering with selected k
  const { labels, clusterColors, actualNumClusters } = useMemo(() => {
    if (scaledData.length === 0 || actualK < 2 || actualK > sampleStats.length) {
      return { labels: [], clusterColors: {}, actualNumClusters: 0 };
    }
    
    const clusterLabels = kMeans(scaledData, actualK, 100, 42);
    
    // Remap cluster labels to be sequential (1, 2, 3, ...) even if some clusters are empty
    const uniqueLabels = Array.from(new Set(clusterLabels)).sort((a, b) => a - b);
    const actualNumClusters = uniqueLabels.length;
    const labelMap: { [key: number]: number } = {};
    uniqueLabels.forEach((oldLabel, idx) => {
      labelMap[oldLabel] = idx + 1; // Map to sequential labels starting from 1
    });
    
    const labels = clusterLabels.map(l => labelMap[l] || 1);
    
    // Color palette for clusters
    const colors = [
      "#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#6A994E", "#7209B7", "#4361EE", "#4CC9F0"
    ];
    const uniqueClusters = Array.from(new Set(labels)).sort((a, b) => a - b);
    const clusterColors: { [key: number]: string } = {};
    uniqueClusters.forEach((cluster, idx) => {
      clusterColors[cluster] = colors[idx % colors.length];
    });
    
    return { labels, clusterColors, actualNumClusters };
  }, [scaledData, actualK, sampleStats.length]);
  
  // Check if sample is current/reference
  const isCurrentSample = (sampleName: string): boolean => {
    return sampleName === 'Ref' || 
           sampleName === 'Current Sample' || 
           sampleName.startsWith('Allel') ||
           sampleName.startsWith('Allele');
  };
  
  const maxCopyNumber = Math.max(...sampleStats.map(s => s.copyNumber), 1);
  const maxLength = Math.max(...sampleStats.map(s => s.length), 1);
  
  if (sampleStats.length < 2) {
    return (
      <div className="info-message">
        Need at least 2 samples for clustering.
      </div>
    );
  }
  
  return (
    <div className="cluster-plot-container">
      <div className="cluster-plot-controls">
        <div className="cluster-control-group">
          <label>Cluster by:</label>
          <div className="cluster-radio-group">
            <label>
              <input
                type="radio"
                value="Copy Number"
                checked={clusterBy === 'Copy Number'}
                onChange={(e) => setClusterBy(e.target.value as 'Copy Number')}
              />
              Copy Number
            </label>
            <label>
              <input
                type="radio"
                value="Length"
                checked={clusterBy === 'Length'}
                onChange={(e) => setClusterBy(e.target.value as 'Length')}
              />
              Length
            </label>
            <label>
              <input
                type="radio"
                value="Both"
                checked={clusterBy === 'Both'}
                onChange={(e) => setClusterBy(e.target.value as 'Both')}
              />
              Both
            </label>
          </div>
        </div>
        
        <div className="cluster-control-group">
          <label htmlFor="num-clusters">
            Number of Clusters (K):
            <span className="cluster-suggested">(Suggested: {suggestedK})</span>
            {actualNumClusters > 0 && actualNumClusters < actualK && (
              <span className="cluster-warning" title={`Only ${actualNumClusters} clusters were created (some clusters were empty)`}>
                ⚠️ Actual: {actualNumClusters}
              </span>
            )}
          </label>
          <input
            id="num-clusters"
            type="number"
            min={2}
            max={Math.min(8, sampleStats.length)}
            value={actualK}
            onChange={(e) => {
              const val = parseInt(e.target.value);
              if (!isNaN(val) && val >= 2 && val <= Math.min(8, sampleStats.length)) {
                setNumClusters(val);
                setUserChangedK(true); // Mark that user manually changed K
              }
            }}
            onBlur={(e) => {
              const val = parseInt(e.target.value);
              if (isNaN(val) || val < 2) {
                setNumClusters(2);
                setUserChangedK(true);
              } else if (val > Math.min(8, sampleStats.length)) {
                setNumClusters(Math.min(8, sampleStats.length));
                setUserChangedK(true);
              }
            }}
            title="Choose number of clusters. Current cluster number is calculated using silhouette score."
          />
        </div>
      </div>
      
      <div className="cluster-plot-scatter">
        <div className="cluster-plot-axis-labels">
          <div className="y-axis-label">
            {clusterBy === 'Both' ? 'Copy Number' : clusterBy}
          </div>
          <div className="cluster-plot-area">
            {sampleStats.map((stat, idx) => {
              const clusterLabel = labels[idx] || 0;
              const color = clusterColors[clusterLabel] || '#9ca3af';
              const isCurrent = isCurrentSample(stat.sample);
              
              let xPercent: number;
              let yPercent: number;
              
              if (clusterBy === 'Both') {
                xPercent = (stat.length / maxLength) * 100;
                yPercent = 100 - (stat.copyNumber / maxCopyNumber) * 100;
              } else if (clusterBy === 'Copy Number') {
                xPercent = (stat.copyNumber / maxCopyNumber) * 100;
                yPercent = 50; // Center vertically for 1D
              } else {
                xPercent = (stat.length / maxLength) * 100;
                yPercent = 50; // Center vertically for 1D
              }
              
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
                  title={`${stat.sample}: Cluster=${clusterLabel}, Copy Number=${stat.copyNumber}, Length=${stat.length}bp`}
                >
                  <span className="cluster-point-label">{stat.sample}</span>
                  {clusterLabel > 0 && (
                    <span className="cluster-number-badge">{clusterLabel}</span>
                  )}
                </div>
              );
            })}
            
            <div className="cluster-grid-lines">
              {[0, 25, 50, 75, 100].map(percent => (
                <React.Fragment key={percent}>
                  <div className="grid-line-vertical" style={{ left: `${percent}%` }} />
                  {clusterBy === 'Both' && (
                    <div className="grid-line-horizontal" style={{ bottom: `${percent}%` }} />
                  )}
                </React.Fragment>
              ))}
            </div>
          </div>
          <div className="x-axis-label">
            {clusterBy === 'Both' ? 'Length (bp)' : clusterBy}
          </div>
        </div>
      </div>
      
      <div className="cluster-plot-legend">
        <div className="legend-clusters">
          <h4>Clusters</h4>
          <div className="cluster-legend-items">
            {Object.entries(clusterColors).map(([cluster, color]) => (
              <div key={cluster} className="cluster-legend-item">
                <span className="cluster-legend-color" style={{ backgroundColor: color }} />
                <span>Cluster {cluster}</span>
              </div>
            ))}
          </div>
        </div>
        <div className="legend-note">
          <span className="legend-current-marker" /> = Current Sample / Reference
        </div>
      </div>
    </div>
  );
});
ClusterPlot.displayName = 'ClusterPlot';

export default CohortAnalysis;

