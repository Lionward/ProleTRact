import React, { useState, useEffect, useRef } from 'react';
import { FilterResponse } from '../types';
import { RegionSearchBar } from './RegionVisualization';
import FileBrowser from './FileBrowser';
import './Sidebar.css';

// Sidebar component with resizable functionality

type AppMode = 'individual' | 'cohort-read' | 'cohort-assembly';

interface SidebarProps {
  mode: AppMode;
  onModeChange: (mode: AppMode) => void;
  vcfPath: string;
  setVcfPath: (path: string) => void;
  onLoadVCF: () => void;
  publicVcfFolder: string;
  setPublicVcfFolder: (path: string) => void;
  onLoadPublicVCF: () => void;
  cohortFolder: string;
  setCohortFolder: (path: string) => void;
  onLoadCohortFolder: () => void;
  availableGenotypes: string[];
  selectedGenotypes: string[];
  onGenotypeChange: (genotypes: string[]) => void;
  loading: boolean;
  // Cohort mode region search
  cohortRegionInput?: string;
  setCohortRegionInput?: (value: string) => void;
  onCohortRegionSubmit?: (selectedRegion?: string) => void;
  cohortAvailableRegions?: string[];
  loadingCohortRegions?: boolean;
  cohortLoadProgress?: { current: number; total: number } | null;
  onOpenSessionManager?: () => void;
}

const Sidebar: React.FC<SidebarProps> = ({
  mode,
  onModeChange,
  vcfPath,
  setVcfPath,
  onLoadVCF,
  publicVcfFolder,
  setPublicVcfFolder,
  onLoadPublicVCF,
  cohortFolder,
  setCohortFolder,
  onLoadCohortFolder,
  availableGenotypes,
  selectedGenotypes,
  onGenotypeChange,
  loading,
  cohortRegionInput = '',
  setCohortRegionInput,
  onCohortRegionSubmit,
  cohortAvailableRegions = [],
  loadingCohortRegions = false,
  cohortLoadProgress = null,
  onOpenSessionManager
}) => {
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [vcfExpanded, setVcfExpanded] = useState(true);
  const [populationExpanded, setPopulationExpanded] = useState(true);
  const [showVcfBrowser, setShowVcfBrowser] = useState(false);
  const [showPublicFolderBrowser, setShowPublicFolderBrowser] = useState(false);
  const [showCohortFolderBrowser, setShowCohortFolderBrowser] = useState(false);
  
  // Resizable sidebar state
  const [sidebarWidth, setSidebarWidth] = useState(() => {
    const saved = localStorage.getItem('sidebar-width');
    return saved ? parseInt(saved, 10) : 320;
  });
  const [isResizing, setIsResizing] = useState(false);
  const sidebarRef = useRef<HTMLDivElement>(null);

  // Auto-collapse sections when files are loaded
  useEffect(() => {
    if (availableGenotypes.length > 0 && vcfPath) {
      setVcfExpanded(false);
    }
  }, [availableGenotypes.length, vcfPath]);

  useEffect(() => {
    if (publicVcfFolder) {
      setPopulationExpanded(false);
    }
  }, [publicVcfFolder]);

  // Save sidebar width to localStorage and update CSS variable
  useEffect(() => {
    if (!sidebarCollapsed) {
      localStorage.setItem('sidebar-width', sidebarWidth.toString());
      document.documentElement.style.setProperty('--sidebar-width', `${sidebarWidth}px`);
    } else {
      document.documentElement.style.setProperty('--sidebar-width', '70px');
    }
  }, [sidebarWidth, sidebarCollapsed]);

  // Initialize CSS variable on mount
  useEffect(() => {
    document.documentElement.style.setProperty('--sidebar-width', sidebarCollapsed ? '70px' : `${sidebarWidth}px`);
  }, []);

  // Handle resize start
  const handleResizeStart = (e: React.MouseEvent) => {
    e.preventDefault();
    setIsResizing(true);
  };

  // Handle resize
  useEffect(() => {
    if (!isResizing) {
      document.body.classList.remove('resizing');
      return;
    }

    document.body.classList.add('resizing');

    const handleMouseMove = (e: MouseEvent) => {
      if (!sidebarRef.current) return;
      
      const newWidth = e.clientX;
      const minWidth = 250;
      const maxWidth = Math.min(window.innerWidth * 0.6, 800);
      
      if (newWidth >= minWidth && newWidth <= maxWidth) {
        setSidebarWidth(newWidth);
      }
    };

    const handleMouseUp = () => {
      setIsResizing(false);
      document.body.classList.remove('resizing');
    };

    document.addEventListener('mousemove', handleMouseMove);
    document.addEventListener('mouseup', handleMouseUp);
    document.body.style.cursor = 'col-resize';
    document.body.style.userSelect = 'none';

    return () => {
      document.removeEventListener('mousemove', handleMouseMove);
      document.removeEventListener('mouseup', handleMouseUp);
      document.body.style.cursor = '';
      document.body.style.userSelect = '';
      document.body.classList.remove('resizing');
    };
  }, [isResizing]);

  return (
    <>
      <div 
        ref={sidebarRef}
        className={`sidebar ${sidebarCollapsed ? 'collapsed' : ''} ${isResizing ? 'resizing' : ''}`}
        style={!sidebarCollapsed ? { width: `${sidebarWidth}px` } : undefined}
      >
        {/* Sidebar Header with Logo */}
        <div className="sidebar-header">
          <div className="sidebar-brand">
            <img 
              src="/ProleTRact_logo.svg" 
              alt="ProleTRact Logo" 
              className={sidebarCollapsed ? "sidebar-logo-collapsed" : "sidebar-logo"}
            />
          </div>
        </div>
        
        {!sidebarCollapsed && (
          <>
            {/* Mode Selector */}
            <div className="sidebar-mode-selector">
              <button
                className={`sidebar-mode-button ${mode === 'individual' ? 'active' : ''}`}
                onClick={() => onModeChange('individual')}
                title="Single VCF file analysis"
              >
                <span className="sidebar-mode-icon">🔬</span>
                <div className="sidebar-mode-content">
                  <span className="sidebar-mode-label">Individual Mode</span>
                  <span className="sidebar-mode-description">Single VCF file</span>
                </div>
              </button>
              <button
                className={`sidebar-mode-button ${mode === 'cohort-read' ? 'active' : ''}`}
                onClick={() => onModeChange('cohort-read')}
                title="Cohort with diploid read-based VCF files"
              >
                <span className="sidebar-mode-icon">👥</span>
                <div className="sidebar-mode-content">
                  <span className="sidebar-mode-label">Cohort (Read-based)</span>
                  <span className="sidebar-mode-description">Diploid VCF files</span>
                </div>
              </button>
              <button
                className={`sidebar-mode-button ${mode === 'cohort-assembly' ? 'active' : ''}`}
                onClick={() => onModeChange('cohort-assembly')}
                title="Cohort with assembly-based single-haplotype VCF files"
              >
                <span className="sidebar-mode-icon">🧬</span>
                <div className="sidebar-mode-content">
                  <span className="sidebar-mode-label">Cohort (Assembly-based)</span>
                  <span className="sidebar-mode-description">Single haplotype</span>
                </div>
              </button>
            </div>

            {mode === 'individual' ? (
              <>
                {/* Individual Mode: VCF File Input */}
                <div className="sidebar-section">
                  <div className="file-accordion">
                    <div
                      className="file-accordion-header modern"
                      onClick={() => setVcfExpanded(!vcfExpanded)}
                    >
                      <div className="accordion-title">
                        <span>Load VCF File</span>
                        {availableGenotypes.length > 0 && (
                          <span className="file-status-badge loaded">✓ Loaded</span>
                        )}
                      </div>
                      <span className="accordion-arrow">{vcfExpanded ? '▼' : '▶'}</span>
                    </div>
                    {vcfExpanded && (
                      <div className="file-accordion-content">
                        <div className="input-group">
                          <label>
                            VCF File Path
                          </label>
                          <div className="input-with-browse">
                            <input
                              type="text"
                              value={vcfPath}
                              onChange={(e) => setVcfPath(e.target.value)}
                              placeholder="/path/to/file.vcf.gz"
                              onKeyPress={(e) => e.key === 'Enter' && onLoadVCF()}
                              className="modern-input"
                            />
                            <button
                              type="button"
                              className="btn-browse"
                              onClick={() => setShowVcfBrowser(true)}
                              title="Browse for VCF file"
                            >
                              Browse
                            </button>
                          </div>
                        </div>
                        <button
                          className="btn btn-primary btn-modern"
                          onClick={() => onLoadVCF()}
                          disabled={loading || !vcfPath}
                        >
                          {loading ? <span className="spinner" /> : null}
                          <span>Load VCF</span>
                        </button>
                      </div>
                    )}
                  </div>
                </div>

                {/* Individual Mode: Public VCF Folder Input (for Population Analysis) - Optional */}
                <div className="sidebar-section">
                  <div className="file-accordion">
                    <div
                      className="file-accordion-header modern"
                      onClick={() => setPopulationExpanded(!populationExpanded)}
                    >
                      <div className="accordion-title">
                        <span>Population Data</span>
                        <span className="label-badge">Optional</span>
                        {publicVcfFolder && (
                          <span className="file-status-badge loaded">✓ Loaded</span>
                        )}
                      </div>
                      <span className="accordion-arrow">{populationExpanded ? '▼' : '▶'}</span>
                    </div>
                    {populationExpanded && (
                      <div className="file-accordion-content">
                        <div className="input-group">
                          <label>
                            Public VCF Folder
                          </label>
                          <div className="input-with-browse">
                            <input
                              type="text"
                              value={publicVcfFolder}
                              onChange={(e) => setPublicVcfFolder(e.target.value)}
                              placeholder="/path/to/population/folder"
                              title="Path to folder containing population VCF files. Sample names will be extracted from VCF headers."
                              className="modern-input"
                            />
                            <button
                              type="button"
                              className="btn-browse"
                              onClick={() => setShowPublicFolderBrowser(true)}
                              title="Browse for folder"
                            >
                              Browse
                            </button>
                          </div>
                          <div className="input-hint">
                            For population comparison analysis
                          </div>
                        </div>
                        <button
                          className="btn btn-secondary btn-modern"
                          onClick={onLoadPublicVCF}
                          disabled={loading || !publicVcfFolder}
                        >
                          {loading ? <span className="spinner" /> : null}
                          <span>Load Population Data</span>
                        </button>
                      </div>
                    )}
                  </div>
                </div>
              </>
            ) : (mode === 'cohort-read' || mode === 'cohort-assembly') ? (
              <>
                {/* Cohort Mode: Cohort Folder Input */}
                <div className="sidebar-section">
                  <div className="file-accordion">
                    <div
                      className="file-accordion-header modern"
                      onClick={() => setVcfExpanded(!vcfExpanded)}
                    >
                      <div className="accordion-title">
                        <span>Load Cohort Folder</span>
                        {cohortFolder && publicVcfFolder && (
                          <span className="file-status-badge loaded">✓ Loaded</span>
                        )}
                      </div>
                      <span className="accordion-arrow">{vcfExpanded ? '▼' : '▶'}</span>
                    </div>
                    {vcfExpanded && (
                      <div className="file-accordion-content">
                        <div className="input-group">
                          <label>
                            Cohort VCF Folder
                          </label>
                          <div className="input-with-browse">
                            <input
                              type="text"
                              value={cohortFolder}
                              onChange={(e) => setCohortFolder(e.target.value)}
                              placeholder="/path/to/cohort/folder"
                              title="Path to folder containing cohort VCF files. Sample names will be extracted from VCF headers."
                              onKeyPress={(e) => e.key === 'Enter' && onLoadCohortFolder()}
                              className="modern-input"
                            />
                            <button
                              type="button"
                              className="btn-browse"
                              onClick={() => setShowCohortFolderBrowser(true)}
                              title="Browse for folder"
                            >
                              Browse
                            </button>
                          </div>
                          <div className="input-hint">
                            Folder should contain multiple VCF files (.vcf.gz or .vcf). Sample names will be extracted from VCF headers.
                          </div>
                        </div>
                        <button
                          className="btn btn-primary btn-modern"
                          onClick={onLoadCohortFolder}
                          disabled={loading || !cohortFolder}
                        >
                          {loading ? <span className="spinner" /> : null}
                          <span>Load Cohort</span>
                        </button>
                        {cohortLoadProgress && cohortLoadProgress.total > 0 && (
                          <div className="cohort-load-progress-container">
                            <div className="cohort-load-progress-bar">
                              <div 
                                className="cohort-load-progress-fill"
                                style={{ 
                                  width: `${(cohortLoadProgress.current / cohortLoadProgress.total) * 100}%` 
                                }}
                              />
                            </div>
                            <div className="cohort-load-progress-text">
                              {cohortLoadProgress.current} / {cohortLoadProgress.total} files loaded
                            </div>
                          </div>
                        )}
                      </div>
                    )}
                  </div>
                </div>

                {/* Cohort Mode: Region Search */}
                {(publicVcfFolder || cohortFolder) && (
                  <div className="sidebar-section">
                    <div className="file-accordion">
                      <div className="file-accordion-header modern" style={{ cursor: 'default' }}>
                        <div className="accordion-title">
                          <span>Enter Region</span>
                        </div>
                      </div>
                      <div className="file-accordion-content">
                        {loadingCohortRegions && (
                          <div style={{ fontSize: '0.9em', color: '#9ca3af', marginBottom: '0.5rem' }}>
                            Loading available regions...
                          </div>
                        )}
                        {!loadingCohortRegions && cohortAvailableRegions.length > 0 && (
                          <div style={{ fontSize: '0.9em', color: '#9ca3af', marginBottom: '0.5rem' }}>
                            {cohortAvailableRegions.length} regions available
                          </div>
                        )}
                        {setCohortRegionInput && onCohortRegionSubmit && (
                          <RegionSearchBar
                            value={cohortRegionInput}
                            onChange={setCohortRegionInput}
                            onSubmit={onCohortRegionSubmit}
                            placeholder="Enter region (e.g., chr1:1000-2000)"
                            availableRegions={cohortAvailableRegions}
                          />
                        )}
                      </div>
                    </div>
                  </div>
                )}
              </>
            ) : null}

            {/* Sessions */}
            {onOpenSessionManager && (
              <div className="sidebar-section">
                <button
                  className="btn btn-secondary btn-modern sidebar-sessions-btn"
                  onClick={onOpenSessionManager}
                  title="Manage saved sessions"
                >
                  <span className="sidebar-sessions-icon">💾</span>
                  <span>Sessions</span>
                </button>
              </div>
            )}

          </>
        )}
        
        {/* Resize Handle */}
        {!sidebarCollapsed && (
          <div
            className="sidebar-resize-handle"
            onMouseDown={handleResizeStart}
            title="Drag to resize sidebar"
          />
        )}
      </div>
      
      {/* Fixed Toggle Buttons - Always visible */}
      <button 
        className={`sidebar-toggle-fixed sidebar-toggle-left ${sidebarCollapsed ? 'hidden' : ''}`}
        onClick={() => setSidebarCollapsed(true)}
        title="Collapse sidebar"
        aria-label="Collapse sidebar"
      >
        ◀
      </button>
      
      <button 
        className={`sidebar-toggle-fixed sidebar-toggle-right ${!sidebarCollapsed ? 'hidden' : ''}`}
        onClick={() => setSidebarCollapsed(false)}
        title="Expand sidebar"
        aria-label="Expand sidebar"
      >
        ▶
      </button>
      
      {/* Overlay for mobile when sidebar is open */}
      {!sidebarCollapsed && <div className="sidebar-overlay" onClick={() => setSidebarCollapsed(true)}></div>}

      {/* File browsers */}
      <FileBrowser
        isOpen={showVcfBrowser}
        onClose={() => setShowVcfBrowser(false)}
        onSelect={(path) => { setVcfPath(path); setShowVcfBrowser(false); }}
        mode="file"
        initialPath={vcfPath ? vcfPath.replace(/\/[^/]+$/, '') : ''}
        title="Select VCF File"
      />
      <FileBrowser
        isOpen={showPublicFolderBrowser}
        onClose={() => setShowPublicFolderBrowser(false)}
        onSelect={(path) => { setPublicVcfFolder(path); setShowPublicFolderBrowser(false); }}
        mode="folder"
        initialPath={publicVcfFolder || ''}
        title="Select Population Folder"
      />
      <FileBrowser
        isOpen={showCohortFolderBrowser}
        onClose={() => setShowCohortFolderBrowser(false)}
        onSelect={(path) => { setCohortFolder(path); setShowCohortFolderBrowser(false); }}
        mode="folder"
        initialPath={cohortFolder || ''}
        title="Select Cohort Folder"
      />
    </>
  );
};

export default Sidebar;

