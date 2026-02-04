import React, { useState, useEffect, useRef, useCallback } from 'react';
import axios from 'axios';
import './App.css';
import Sidebar from './components/Sidebar';
import MainContent from './components/MainContent';
import FloatingNavigation from './components/FloatingNavigation';
import SessionManager from './components/SessionManager';
import FilterPanel from './components/FilterPanel';
import { useSession } from './contexts/SessionContext';
import { VCFData, FilterResponse } from './types';

// API base URL - works both locally and through SSH forwarding
// When accessing via SSH forwarding, use localhost since ports are forwarded
const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

type AppMode = 'individual' | 'cohort-read' | 'cohort-assembly';

// Cache structure for mode-specific state
interface ModeCache {
  // Individual mode cache
  individual: {
    vcfPath: string;
    vcfData: VCFData | null;
    selectedGenotypes: string[];
    availableGenotypes: string[];
    availableChromosomes: string[];
    regions: FilterResponse | null;
    currentPage: number;
    selectedRegion: string;
    publicVcfFolder: string;
  };
  // Cohort mode cache (shared between read and assembly)
  'cohort-read': {
    cohortFolder: string;
    publicVcfFolder: string;
    selectedRegion: string;
    cohortRegionInput: string;
    cohortAvailableRegions: string[];
  };
  'cohort-assembly': {
    cohortFolder: string;
    publicVcfFolder: string;
    selectedRegion: string;
    cohortRegionInput: string;
    cohortAvailableRegions: string[];
  };
}

const App: React.FC = () => {
  const [mode, setMode] = useState<AppMode>('individual');
  const [vcfPath, setVcfPath] = useState<string>('');
  const [vcfData, setVcfData] = useState<VCFData | null>(null);
  const [selectedGenotypes, setSelectedGenotypes] = useState<string[]>([]);
  const [availableGenotypes, setAvailableGenotypes] = useState<string[]>([]);
  const [availableChromosomes, setAvailableChromosomes] = useState<string[]>([]);
  const [regions, setRegions] = useState<FilterResponse | null>(null);
  const [currentPage, setCurrentPage] = useState(0);
  const [loading, setLoading] = useState(false);
  const [selectedRegion, setSelectedRegion] = useState<string>('');
  const [publicVcfFolder, setPublicVcfFolder] = useState<string>('');
  const [cohortFolder, setCohortFolder] = useState<string>('');
  // Cohort mode region search state
  const [cohortRegionInput, setCohortRegionInput] = useState<string>('');
  const [cohortAvailableRegions, setCohortAvailableRegions] = useState<string[]>([]);
  const [loadingCohortRegions, setLoadingCohortRegions] = useState(false);
  const [cohortLoadProgress, setCohortLoadProgress] = useState<{ current: number; total: number } | null>(null);
  const [showSessionManager, setShowSessionManager] = useState(false);
  const [showFilterPanel, setShowFilterPanel] = useState(false);
  const [filterClearKey, setFilterClearKey] = useState(0);

  const { currentSession, updateCurrentSession } = useSession();
  const isRestoringFromSession = useRef(false);

  // Get current app state for session save (used when no currentSession exists)
  const getSessionData = useCallback(() => ({
    vcfPath: vcfPath || '',
    selectedRegion: selectedRegion || '',
    selectedGenotypes: selectedGenotypes || [],
    publicVcfFolder: publicVcfFolder || '',
    cohortFolder: cohortFolder || '',
    mode: mode as 'individual' | 'cohort-read' | 'cohort-assembly',
  }), [vcfPath, selectedRegion, selectedGenotypes, publicVcfFolder, cohortFolder, mode]);

  // Restore app state when user loads a session from Session Manager
  const handleSessionLoadRequested = useCallback((session: { vcfPath?: string; selectedRegion?: string; selectedGenotypes?: string[]; publicVcfFolder?: string; cohortFolder?: string; mode?: string }) => {
    if (isRestoringFromSession.current) return;
    isRestoringFromSession.current = true;
    const targetMode = (session.mode === 'cohort-read' || session.mode === 'cohort-assembly') ? session.mode as AppMode : 'individual';
    setMode(targetMode);
    setVcfPath(session.vcfPath || '');
    setSelectedRegion(session.selectedRegion || '');
    setSelectedGenotypes(session.selectedGenotypes || []);
    setPublicVcfFolder(session.publicVcfFolder || '');
    setCohortFolder(session.cohortFolder || '');
    setCohortRegionInput(session.selectedRegion || '');
    const folder = session.publicVcfFolder || session.cohortFolder || '';
    if (session.vcfPath && targetMode === 'individual') {
      loadVCF(session.vcfPath);
    } else if (folder && (targetMode === 'cohort-read' || targetMode === 'cohort-assembly')) {
      loadCohortFolder(folder);
    }
    if (session.selectedRegion && folder && (targetMode === 'cohort-read' || targetMode === 'cohort-assembly')) {
      fetchCohortRegions(folder);
    }
    isRestoringFromSession.current = false;
  }, []);

  // Sync app state to session context when user makes changes (for save to capture current state)
  useEffect(() => {
    if (isRestoringFromSession.current || !currentSession) return;
    updateCurrentSession(getSessionData());
  }, [mode, selectedRegion, vcfPath, selectedGenotypes, publicVcfFolder, cohortFolder, currentSession?.id]);

  // Cache for mode-specific state
  const modeCache = useRef<ModeCache>({
    individual: {
      vcfPath: '',
      vcfData: null,
      selectedGenotypes: [],
      availableGenotypes: [],
      availableChromosomes: [],
      regions: null,
      currentPage: 0,
      selectedRegion: '',
      publicVcfFolder: ''
    },
    'cohort-read': {
      cohortFolder: '',
      publicVcfFolder: '',
      selectedRegion: '',
      cohortRegionInput: '',
      cohortAvailableRegions: []
    },
    'cohort-assembly': {
      cohortFolder: '',
      publicVcfFolder: '',
      selectedRegion: '',
      cohortRegionInput: '',
      cohortAvailableRegions: []
    }
  });

  const loadVCF = async (overridePath?: string) => {
    const pathToUse = (typeof overridePath === 'string' ? overridePath : null) ?? vcfPath;
    // #region agent log
    fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:95',message:'loadVCF called',data:{vcfPath:pathToUse},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'VCF'})}).catch(()=>{});
    // #endregion
    if (!pathToUse) {
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:98',message:'loadVCF early return: no vcfPath',data:{},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'VCF'})}).catch(()=>{});
      // #endregion
      return;
    }
    
    setLoading(true);
    try {
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:103',message:'Calling API to load VCF',data:{vcfPath,apiUrl:`${API_BASE}/api/vcf/load`},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'VCF'})}).catch(()=>{});
      // #endregion
      const response = await axios.post(`${API_BASE}/api/vcf/load`, {
        vcf_path: pathToUse
      });
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:108',message:'VCF load response received',data:{totalRegions:response.data.total_regions,availableGenotypes:response.data.available_genotypes},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'VCF'})}).catch(()=>{});
      // #endregion
      
      setVcfData({
        path: pathToUse,
        totalRegions: response.data.total_regions ?? null,
        availableGenotypes: response.data.available_genotypes ?? [],
        availableChromosomes: response.data.available_chromosomes ?? []
      });
      
      setAvailableGenotypes(response.data.available_genotypes ?? []);
      setAvailableChromosomes(response.data.available_chromosomes ?? []);
      setSelectedGenotypes(response.data.available_genotypes ?? []);
      
      // Load first page (builds index on first use when lazy-loaded)
      await filterRegions(response.data.available_genotypes ?? [], 0, undefined, false, pathToUse);
    } catch (error: any) {
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:124',message:'Error loading VCF',data:{errorMessage:error.message,errorDetail:error.response?.data?.detail,statusCode:error.response?.status},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'VCF'})}).catch(()=>{});
      // #endregion
      console.error('Error loading VCF:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const filterRegions = async (genotypes: string[], page: number, preserveSelectedRegion?: string, selectLastRegion?: boolean, overrideVcfPath?: string) => {
    const pathToUse = overrideVcfPath ?? vcfPath;
    if (!pathToUse) return;
    
    setLoading(true);
    try {
      const response = await axios.post<FilterResponse>(`${API_BASE}/api/vcf/filter`, {
        vcf_path: pathToUse,
        genotype_filter: genotypes.length === availableGenotypes.length ? null : genotypes,
        page: page,
        page_size: 50
      });
      
      setRegions(response.data);
      setCurrentPage(page);
      // Lazy load: first filter builds index and returns available_genotypes
      if (response.data.available_genotypes?.length) {
        setAvailableGenotypes(response.data.available_genotypes);
        setSelectedGenotypes(response.data.available_genotypes);
      }
      
      // If we're preserving a selected region and it's on this page, keep it
      // If selectLastRegion is true, select the last region on the page
      // Otherwise, select the first region on the page
      if (preserveSelectedRegion) {
        const regionOnPage = response.data.records.find(r => r.region === preserveSelectedRegion);
        if (regionOnPage) {
          setSelectedRegion(preserveSelectedRegion);
        } else if (response.data.records.length > 0) {
          setSelectedRegion(response.data.records[0].region);
        }
      } else if (selectLastRegion && response.data.records.length > 0) {
        setSelectedRegion(response.data.records[response.data.records.length - 1].region);
      } else if (response.data.records.length > 0) {
        setSelectedRegion(response.data.records[0].region);
      }
    } catch (error: any) {
      console.error('Error filtering regions:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handleGenotypeChange = (genotypes: string[]) => {
    setSelectedGenotypes(genotypes);
    filterRegions(genotypes, 0);
  };

  // Advanced filter function
  const handleAdvancedFilter = async (criteria: any, page: number = 0) => {
    // #region agent log
    fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:167',message:'handleAdvancedFilter called',data:{criteria,page,vcfPath},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'C'})}).catch(()=>{});
    // #endregion
    if (!vcfPath) {
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:170',message:'handleAdvancedFilter early return: no vcfPath',data:{},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'C'})}).catch(()=>{});
      // #endregion
      return;
    }
    
    setLoading(true);
    try {
      // Collect annotations from localStorage for filtering
      let annotatedRegions = new Set<string>();
      if (criteria.hasAnnotations || criteria.annotationTags) {
        try {
          // Try both storage keys (for compatibility)
          const stored1 = localStorage.getItem('proletract_annotations');
          if (stored1) {
            const allAnnotations = JSON.parse(stored1);
            if (Array.isArray(allAnnotations)) {
              allAnnotations.forEach((ann: any) => {
                if (ann.region) {
                  if (criteria.hasAnnotations) {
                    annotatedRegions.add(ann.region);
                  } else if (criteria.annotationTags && criteria.annotationTags.length > 0) {
                    const annTags = ann.tags || [];
                    if (criteria.annotationTags.some((tag: string) => annTags.includes(tag))) {
                      annotatedRegions.add(ann.region);
                    }
                  }
                }
              });
            }
          }
          // Also check individual region annotations
          for (let i = 0; i < localStorage.length; i++) {
            const key = localStorage.key(i);
            if (key && key.startsWith('annotation_')) {
              const ann = JSON.parse(localStorage.getItem(key) || '{}');
              if (ann.region) {
                if (criteria.hasAnnotations) {
                  annotatedRegions.add(ann.region);
                } else if (criteria.annotationTags && criteria.annotationTags.length > 0) {
                  const annTags = ann.tags || [];
                  if (criteria.annotationTags.some((tag: string) => annTags.includes(tag))) {
                    annotatedRegions.add(ann.region);
                  }
                }
              }
            }
          }
        } catch (err) {
          console.error('Error loading annotations from localStorage:', err);
        }
      }
      
      const apiPayload = {
        vcf_path: vcfPath,
        motif_size_min: criteria.motifSizeMin,
        motif_size_max: criteria.motifSizeMax,
        cn_min: criteria.cnMin,
        cn_max: criteria.cnMax,
        chromosomes: criteria.chromosomes,
        genotypes: criteria.genotypes,
        pathogenic_only: criteria.pathogenicOnly || false,
        has_annotations: false, // We'll filter on frontend for annotations
        annotation_tags: criteria.annotationTags,
        annotated_regions: Array.from(annotatedRegions), // Send annotated regions to backend
        page: page,
        page_size: 50
      };
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:194',message:'Making API call to filter-advanced',data:{apiPayload,url:`${API_BASE}/api/vcf/filter-advanced`},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      
      const response = await axios.post<FilterResponse>(`${API_BASE}/api/vcf/filter-advanced`, apiPayload);
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:199',message:'API response received',data:{totalMatching:response.data.total_matching,recordsCount:response.data.records.length,totalPages:response.data.total_pages},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      
      setRegions(response.data);
      setCurrentPage(page);
      if (response.data.available_genotypes?.length) {
        setAvailableGenotypes(response.data.available_genotypes);
        setSelectedGenotypes(response.data.available_genotypes);
      }
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:266',message:'State updated with filtered regions',data:{totalMatching:response.data.total_matching,recordsCount:response.data.records.length,hasRecords:response.data.records.length>0},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
      // #endregion
      
      // Update selected genotypes if filter includes genotypes
      if (criteria.genotypes && criteria.genotypes.length > 0) {
        setSelectedGenotypes(criteria.genotypes);
      }
      
      // Select first region on page if available
      if (response.data.records.length > 0) {
        setSelectedRegion(response.data.records[0].region);
        // #region agent log
        fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:279',message:'Selected first region',data:{selectedRegion:response.data.records[0].region},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
      } else {
        // No regions match the filter - clear selection and show message
        setSelectedRegion('');
        if (response.data.total_matching === 0) {
          // Show informative message about why no results
          const message = criteria.pathogenicOnly 
            ? 'No pathogenic regions found. All regions have copy numbers below the pathogenic threshold.'
            : 'No regions match the selected filter criteria.';
          alert(message);
        }
        // #region agent log
        fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:287',message:'No regions match filter',data:{totalMatching:response.data.total_matching,pathogenicOnly:criteria.pathogenicOnly},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'E'})}).catch(()=>{});
        // #endregion
      }
    } catch (error: any) {
      // #region agent log
      fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'App.tsx:219',message:'Error in handleAdvancedFilter',data:{errorMessage:error.message,errorDetail:error.response?.data?.detail,statusCode:error.response?.status},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'D'})}).catch(()=>{});
      // #endregion
      console.error('Error applying advanced filter:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const handlePageChange = (newPage: number) => {
    filterRegions(selectedGenotypes, newPage);
  };

  // Jump to a specific region number (1-indexed)
  const handleJumpToRegion = async (regionNumber: number) => {
    if (!regions || !vcfPath) return;
    
    // Convert to 0-indexed
    const regionIndex = regionNumber - 1;
    
    // Validate range
    if (regionIndex < 0 || regionIndex >= regions.total_matching) {
      alert(`Region number must be between 1 and ${regions.total_matching.toLocaleString()}`);
      return;
    }
    
    try {
      const genotypeFilter = selectedGenotypes.length === availableGenotypes.length 
        ? null 
        : selectedGenotypes.join(',');
      
      const response = await axios.get(`${API_BASE}/api/vcf/region-by-index`, {
        params: {
          vcf_path: vcfPath,
          region_index: regionIndex,
          genotype_filter: genotypeFilter,
          page_size: 50
        }
      });
      
      if (response.data.success && response.data.region) {
        const targetPage = response.data.page;
        // Navigate to that page and select the region
        await filterRegions(selectedGenotypes, targetPage, response.data.region);
      }
    } catch (error: any) {
      console.error('Error finding region by index:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    }
  };

  // Navigate to next region (single step)
  const handleNextRegion = () => {
    if (!regions || !selectedRegion) return;
    
    const currentIndex = regions.records.findIndex(r => r.region === selectedRegion);
    
    if (currentIndex === -1) {
      if (regions.records.length > 0) {
        setSelectedRegion(regions.records[0].region);
      }
      return;
    }
    
    if (currentIndex < regions.records.length - 1) {
      setSelectedRegion(regions.records[currentIndex + 1].region);
    } else if (currentPage < regions.total_pages - 1) {
      filterRegions(selectedGenotypes, currentPage + 1);
    }
  };

  // Navigate to previous region (single step)
  const handlePreviousRegion = async () => {
    if (!regions || !selectedRegion) return;
    
    const currentIndex = regions.records.findIndex(r => r.region === selectedRegion);
    
    if (currentIndex === -1) {
      if (regions.records.length > 0) {
        setSelectedRegion(regions.records[regions.records.length - 1].region);
      }
      return;
    }
    
    if (currentIndex > 0) {
      setSelectedRegion(regions.records[currentIndex - 1].region);
    } else if (currentPage > 0) {
      await filterRegions(selectedGenotypes, currentPage - 1, undefined, true);
    }
  };

  // Find and navigate to the page containing a specific region
  const findRegionPage = async (region: string) => {
    if (!vcfPath || !regions) return;
    
    // First check if the region is in the current page
    const currentRegionIndex = regions.records.findIndex(r => r.region === region);
    if (currentRegionIndex !== -1) {
      // Region is already on the current page, no need to change
      return;
    }
    
    // Region is not on current page, find which page it's on
    try {
      const genotypeFilter = selectedGenotypes.length === availableGenotypes.length 
        ? null 
        : selectedGenotypes.join(',');
      
      const response = await axios.get(`${API_BASE}/api/vcf/region-page`, {
        params: {
          vcf_path: vcfPath,
          region: region,
          genotype_filter: genotypeFilter,
          page_size: 50
        }
      });
      
      if (response.data.success) {
        const targetPage = response.data.page;
        // Navigate to that page, preserving the selected region
        await filterRegions(selectedGenotypes, targetPage, region);
      }
    } catch (error: any) {
      console.error('Error finding region page:', error);
      // Silently fail - region might not be in filtered results
    }
  };

  // Handle region selection - find its page and navigate to it
  const handleRegionSelect = async (region: string) => {
    setSelectedRegion(region);
    // Update cohort region input if in cohort mode
    if (isCohortMode) {
      setCohortRegionInput(region);
    }
    if (mode === 'individual' && vcfPath && regions) {
      await findRegionPage(region);
    }
  };

  const loadPublicVCF = async () => {
    if (!publicVcfFolder) return;
    
    setLoading(true);
    try {
      const response = await axios.post(`${API_BASE}/api/population/load`, {
        folder_path: publicVcfFolder
      });
      
      if (response.data.success) {
        alert(`Loaded ${response.data.file_count} population VCF files`);
      }
    } catch (error: any) {
      console.error('Error loading population VCF files:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const loadCohortFolder = async (overrideFolder?: string) => {
    const folderToUse = overrideFolder ?? cohortFolder;
    if (!folderToUse) return;
    
    // Skip redundant load when same folder is already loaded (prevents error on second press)
    const norm = (p: string) => (p || '').replace(/\\/g, '/').replace(/\/+$/, '');
    if (publicVcfFolder && norm(publicVcfFolder) === norm(folderToUse)) {
      alert('Cohort folder already loaded. Enter a different path to load another cohort.');
      return;
    }
    
    setLoading(true);
    setCohortLoadProgress({ current: 0, total: 1 }); // Show initial progress
    try {
      const response = await axios.post(`${API_BASE}/api/population/load`, {
        folder_path: folderToUse
      });
      
      if (response.data.success) {
        // Store the cohort folder for use in cohort mode
        setPublicVcfFolder(folderToUse);
        setCohortFolder(folderToUse);
        // Update progress to show completion
        setCohortLoadProgress({ current: response.data.file_count, total: response.data.file_count });
        alert(`Loaded ${response.data.file_count} cohort VCF files`);
        // Fetch available regions for the loaded cohort
        fetchCohortRegions(folderToUse);
        // Clear progress after a short delay
        setTimeout(() => setCohortLoadProgress(null), 3000);
      }
    } catch (error: any) {
      console.error('Error loading cohort VCF files:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
      setCohortLoadProgress(null);
    } finally {
      setLoading(false);
    }
  };

  // Fetch available regions for cohort mode
  const fetchCohortRegions = async (folderPath: string) => {
    if (!folderPath) {
      setCohortAvailableRegions([]);
      return;
    }

    setLoadingCohortRegions(true);
    try {
      const response = await axios.get(`${API_BASE}/api/population/regions`, {
        params: { folder_path: folderPath }
      });
      const regions = response.data.regions || [];
      setCohortAvailableRegions(regions);
      
      // Auto-select the first region if available and no region is selected
      if (regions.length > 0 && !selectedRegion) {
        const firstRegion = regions[0];
        setSelectedRegion(firstRegion);
        setCohortRegionInput(firstRegion);
      }
    } catch (err: any) {
      console.error('Failed to fetch regions for autocomplete:', err);
      setCohortAvailableRegions([]);
    } finally {
      setLoadingCohortRegions(false);
    }
  };

  // Handle cohort region submit
  const handleCohortRegionSubmit = (selectedRegionFromAutocomplete?: string) => {
    const cohortFolderPath = publicVcfFolder || cohortFolder;
    if (!cohortFolderPath) {
      alert('Please load a cohort folder first');
      return;
    }
    
    const regionToUse = selectedRegionFromAutocomplete || cohortRegionInput.trim();
    
    // Validate region format: chr:start-end
    const regionPattern = /^([\w]+):(\d+)-(\d+)$/;
    const match = regionToUse.match(regionPattern);
    if (match) {
      setSelectedRegion(regionToUse);
      setCohortRegionInput(regionToUse);
    } else {
      alert('Please enter a valid region format: chr:start-end (e.g., chr1:1000-2000)');
    }
  };

  // Helper to check if current mode is any cohort mode
  const isCohortMode = mode === 'cohort-read' || mode === 'cohort-assembly';

  // Keep individual mode cache updated
  useEffect(() => {
    if (mode === 'individual') {
      modeCache.current.individual = {
        vcfPath,
        vcfData,
        selectedGenotypes,
        availableGenotypes,
        availableChromosomes,
        regions,
        currentPage,
        selectedRegion,
        publicVcfFolder
      };
    }
  }, [mode, vcfPath, vcfData, selectedGenotypes, availableGenotypes, availableChromosomes, regions, currentPage, selectedRegion, publicVcfFolder]);

  // Keep cohort mode cache updated
  useEffect(() => {
    if (isCohortMode) {
      modeCache.current[mode] = {
        cohortFolder,
        publicVcfFolder,
        selectedRegion,
        cohortRegionInput,
        cohortAvailableRegions
      };
    }
  }, [mode, cohortFolder, publicVcfFolder, selectedRegion, cohortRegionInput, cohortAvailableRegions, isCohortMode]);

  // Fetch regions when cohort folder changes (only if not already cached)
  useEffect(() => {
    if (isCohortMode) {
      const cohortFolderPath = publicVcfFolder || cohortFolder;
      if (cohortFolderPath && cohortAvailableRegions.length === 0) {
        fetchCohortRegions(cohortFolderPath);
      } else if (!cohortFolderPath) {
        setCohortAvailableRegions([]);
        setCohortRegionInput('');
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [mode, publicVcfFolder, cohortFolder, isCohortMode]);

  // Save current mode state to cache before switching
  const saveCurrentModeState = () => {
    if (mode === 'individual') {
      modeCache.current.individual = {
        vcfPath,
        vcfData,
        selectedGenotypes,
        availableGenotypes,
        availableChromosomes,
        regions,
        currentPage,
        selectedRegion,
        publicVcfFolder
      };
    } else if (isCohortMode) {
      modeCache.current[mode] = {
        cohortFolder,
        publicVcfFolder,
        selectedRegion,
        cohortRegionInput,
        cohortAvailableRegions
      };
    }
  };

  // Restore state from cache when switching to a mode
  const restoreModeState = (targetMode: AppMode) => {
    if (targetMode === 'individual') {
      const cached = modeCache.current.individual;
      setVcfPath(cached.vcfPath);
      setVcfData(cached.vcfData);
      setSelectedGenotypes(cached.selectedGenotypes);
      setAvailableGenotypes(cached.availableGenotypes);
      setAvailableChromosomes(cached.availableChromosomes);
      setRegions(cached.regions);
      setCurrentPage(cached.currentPage);
      setSelectedRegion(cached.selectedRegion);
      setPublicVcfFolder(cached.publicVcfFolder);
    } else if (targetMode === 'cohort-read' || targetMode === 'cohort-assembly') {
      const cached = modeCache.current[targetMode];
      setCohortFolder(cached.cohortFolder);
      setPublicVcfFolder(cached.publicVcfFolder);
      setSelectedRegion(cached.selectedRegion);
      setCohortRegionInput(cached.cohortRegionInput);
      setCohortAvailableRegions(cached.cohortAvailableRegions);
      
      // If cohort folder is loaded, fetch regions if not already cached
      const cohortFolderPath = cached.publicVcfFolder || cached.cohortFolder;
      if (cohortFolderPath && cached.cohortAvailableRegions.length === 0) {
        fetchCohortRegions(cohortFolderPath);
      }
    }
  };

  // Reset state when switching modes
  const handleModeChange = (newMode: AppMode) => {
    if (newMode !== mode) {
      // Save current mode state before switching
      saveCurrentModeState();
      
      // Switch mode
      setMode(newMode);
      
      // Restore the new mode's cached state
      restoreModeState(newMode);
    }
  };

  return (
    <div className="app-container">
      {showSessionManager && (
        <SessionManager
          onClose={() => setShowSessionManager(false)}
          onSessionLoadRequested={handleSessionLoadRequested}
          getSessionData={getSessionData}
        />
      )}
      <Sidebar
        mode={mode}
        onModeChange={handleModeChange}
        vcfPath={vcfPath}
        setVcfPath={setVcfPath}
        onLoadVCF={loadVCF}
        publicVcfFolder={publicVcfFolder}
        setPublicVcfFolder={setPublicVcfFolder}
        onLoadPublicVCF={loadPublicVCF}
        cohortFolder={cohortFolder}
        setCohortFolder={setCohortFolder}
        onLoadCohortFolder={loadCohortFolder}
        availableGenotypes={availableGenotypes}
        selectedGenotypes={selectedGenotypes}
        onGenotypeChange={handleGenotypeChange}
        loading={loading}
        cohortRegionInput={cohortRegionInput}
        setCohortRegionInput={setCohortRegionInput}
        onCohortRegionSubmit={handleCohortRegionSubmit}
        cohortAvailableRegions={cohortAvailableRegions}
        loadingCohortRegions={loadingCohortRegions}
        cohortLoadProgress={cohortLoadProgress}
        onOpenSessionManager={() => setShowSessionManager(true)}
      />
      <div className="app-wrapper">
        <MainContent
          mode={mode}
          selectedRegion={selectedRegion}
          vcfPath={vcfPath}
          loading={loading}
          publicVcfFolder={publicVcfFolder}
          cohortFolder={cohortFolder}
          regions={mode === 'individual' ? regions : undefined}
          onRegionSelect={handleRegionSelect}
          onOpenSessionManager={() => setShowSessionManager(true)}
          getSessionData={getSessionData}
          onClearFilter={mode === 'individual' && vcfPath ? async () => {
            localStorage.removeItem('filter_criteria_current');
            setFilterClearKey(k => k + 1);
            await handleAdvancedFilter({}, 0);
          } : undefined}
        />
        {mode === 'individual' && (
          <>
            {/* Filter button - right side, vertically centered */}
            {vcfPath && availableGenotypes.length > 0 && (
              <button
                className="filter-fab"
                onClick={() => setShowFilterPanel(true)}
                title="Open filters"
                aria-label="Open filters"
              >
                <svg className="filter-icon" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                  <polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3" />
                </svg>
              </button>
            )}
            <FilterPanel
              key={filterClearKey}
              isOpen={showFilterPanel}
              onClose={() => setShowFilterPanel(false)}
              vcfPath={vcfPath}
              availableGenotypes={availableGenotypes}
              selectedGenotypes={selectedGenotypes}
              availableChromosomes={availableChromosomes}
              loading={loading}
              onApplyFilter={(criteria) => handleAdvancedFilter(criteria, 0)}
            />
            <FloatingNavigation
              regions={regions}
              currentPage={currentPage}
              onPageChange={handlePageChange}
              selectedRegion={selectedRegion}
              onNextRegion={handleNextRegion}
              onPreviousRegion={handlePreviousRegion}
              onJumpToRegion={handleJumpToRegion}
              loading={loading}
            />
          </>
        )}
      </div>
    </div>
  );
};

export default App;
