import React, { useState, useEffect, useRef } from 'react';
import axios from 'axios';
import './App.css';
import Sidebar from './components/Sidebar';
import MainContent from './components/MainContent';
import FloatingNavigation from './components/FloatingNavigation';
import SessionManager from './components/SessionManager';
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
  const [pathogenicRegions, setPathogenicRegions] = useState<Set<string>>(new Set());
  const [pathogenicFilterActive, setPathogenicFilterActive] = useState(false);

  // Cache for mode-specific state
  const modeCache = useRef<ModeCache>({
    individual: {
      vcfPath: '',
      vcfData: null,
      selectedGenotypes: [],
      availableGenotypes: [],
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

  const loadVCF = async () => {
    if (!vcfPath) return;
    
    setLoading(true);
    try {
      const response = await axios.post(`${API_BASE}/api/vcf/load`, {
        vcf_path: vcfPath
      });
      
      setVcfData({
        path: vcfPath,
        totalRegions: response.data.total_regions,
        availableGenotypes: response.data.available_genotypes
      });
      
      setAvailableGenotypes(response.data.available_genotypes);
      setSelectedGenotypes(response.data.available_genotypes);
      
      // Load first page
      await filterRegions(response.data.available_genotypes, 0);
    } catch (error: any) {
      console.error('Error loading VCF:', error);
      alert(`Error: ${error.response?.data?.detail || error.message}`);
    } finally {
      setLoading(false);
    }
  };

  const filterRegions = async (genotypes: string[], page: number, preserveSelectedRegion?: string, selectLastRegion?: boolean) => {
    if (!vcfPath) return;
    
    setLoading(true);
    try {
      const response = await axios.post<FilterResponse>(`${API_BASE}/api/vcf/filter`, {
        vcf_path: vcfPath,
        genotype_filter: genotypes.length === availableGenotypes.length ? null : genotypes,
        page: page,
        page_size: 50
      });
      
      setRegions(response.data);
      setCurrentPage(page);
      
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
    
    // Filter regions if pathogenic filter is active
    const filteredRecords = pathogenicFilterActive && pathogenicRegions.size > 0
      ? regions.records.filter(r => pathogenicRegions.has(r.region))
      : regions.records;
    
    const currentIndex = filteredRecords.findIndex(r => r.region === selectedRegion);
    
    if (currentIndex === -1) {
      // Selected region not found, go to first region of filtered list
      if (filteredRecords.length > 0) {
        setSelectedRegion(filteredRecords[0].region);
      }
      return;
    }
    
    // If not at the end of filtered list, move to next region
    if (currentIndex < filteredRecords.length - 1) {
      setSelectedRegion(filteredRecords[currentIndex + 1].region);
    } else {
      // At end of filtered list, try to find next region in next pages
      // For now, just stay at last region if filter is active
      if (!pathogenicFilterActive && currentPage < regions.total_pages - 1) {
        filterRegions(selectedGenotypes, currentPage + 1);
      }
    }
  };

  // Navigate to previous region (single step)
  const handlePreviousRegion = async () => {
    if (!regions || !selectedRegion) return;
    
    // Filter regions if pathogenic filter is active
    const filteredRecords = pathogenicFilterActive && pathogenicRegions.size > 0
      ? regions.records.filter(r => pathogenicRegions.has(r.region))
      : regions.records;
    
    const currentIndex = filteredRecords.findIndex(r => r.region === selectedRegion);
    
    if (currentIndex === -1) {
      // Selected region not found, go to last region of filtered list
      if (filteredRecords.length > 0) {
        setSelectedRegion(filteredRecords[filteredRecords.length - 1].region);
      }
      return;
    }
    
    // If not at the beginning of filtered list, move to previous region
    if (currentIndex > 0) {
      setSelectedRegion(filteredRecords[currentIndex - 1].region);
    } else {
      // At beginning of filtered list, try to find previous region in previous pages
      // For now, just stay at first region if filter is active
      if (!pathogenicFilterActive && currentPage > 0) {
        await filterRegions(selectedGenotypes, currentPage - 1, undefined, true);
      }
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

  const loadCohortFolder = async () => {
    if (!cohortFolder) return;
    
    setLoading(true);
    setCohortLoadProgress({ current: 0, total: 1 }); // Show initial progress
    try {
      const response = await axios.post(`${API_BASE}/api/population/load`, {
        folder_path: cohortFolder
      });
      
      if (response.data.success) {
        // Store the cohort folder for use in cohort mode
        setPublicVcfFolder(cohortFolder);
        // Update progress to show completion
        setCohortLoadProgress({ current: response.data.file_count, total: response.data.file_count });
        alert(`Loaded ${response.data.file_count} cohort VCF files`);
        // Fetch available regions for the loaded cohort
        fetchCohortRegions(cohortFolder);
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
        regions,
        currentPage,
        selectedRegion,
        publicVcfFolder
      };
    }
  }, [mode, vcfPath, vcfData, selectedGenotypes, availableGenotypes, regions, currentPage, selectedRegion, publicVcfFolder]);

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
        <SessionManager onClose={() => setShowSessionManager(false)} />
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
      />
      <div className="app-wrapper">
        <MainContent
          mode={mode}
          selectedRegion={selectedRegion}
          vcfPath={vcfPath}
          loading={loading}
          publicVcfFolder={publicVcfFolder}
          cohortFolder={cohortFolder}
          onRegionSelect={handleRegionSelect}
          onOpenSessionManager={() => setShowSessionManager(true)}
          onPathogenicRegionsChange={(regions, filterActive) => {
            setPathogenicRegions(regions);
            setPathogenicFilterActive(filterActive);
          }}
        />
        {mode === 'individual' && (
          <FloatingNavigation
            regions={pathogenicFilterActive && pathogenicRegions.size > 0 && regions ? (() => {
              const filteredRecords = regions.records.filter(r => pathogenicRegions.has(r.region));
              const totalPages = Math.max(1, Math.ceil(pathogenicRegions.size / 50));
              return {
                ...regions,
                records: filteredRecords,
                total_matching: pathogenicRegions.size,
                total_pages: totalPages,
                current_page: currentPage
              };
            })() : regions}
            currentPage={currentPage}
            onPageChange={handlePageChange}
            selectedRegion={selectedRegion}
            onNextRegion={handleNextRegion}
            onPreviousRegion={handlePreviousRegion}
            onJumpToRegion={handleJumpToRegion}
            loading={loading}
          />
        )}
      </div>
    </div>
  );
};

export default App;
