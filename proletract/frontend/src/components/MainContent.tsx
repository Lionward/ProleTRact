import React, { useState, useRef } from 'react';
import './MainContent.css';
import RegionVisualization from './RegionVisualization';
import StatisticsDashboard from './StatisticsDashboard';
import CohortModeView from './CohortModeView';
import BookmarksPanel from './BookmarksPanel';
import KeyboardShortcutsHelp from './KeyboardShortcutsHelp';
import Tutorial from './Tutorial';
import InstructionsPage from './InstructionsPage';
import FeaturesPage from './FeaturesPage';
import TandemRepeatAnimation from './TandemRepeatAnimation';
import { useKeyboardShortcuts, createAppShortcuts } from '../hooks/useKeyboardShortcuts';
import { useTheme } from '../contexts/ThemeContext';
import { FilterResponse } from '../types';

type AppMode = 'individual' | 'cohort-read' | 'cohort-assembly';

interface MainContentProps {
  mode: AppMode;
  selectedRegion: string;
  vcfPath: string;
  loading: boolean;
  publicVcfFolder?: string;
  cohortFolder?: string;
  regions?: FilterResponse | null;
  onRegionSelect?: (region: string) => void;
  onOpenSessionManager?: () => void;
  onClearFilter?: () => void;
  getSessionData?: () => { vcfPath: string; selectedRegion: string; selectedGenotypes: string[]; publicVcfFolder: string; cohortFolder: string; mode: string };
}

const MainContent: React.FC<MainContentProps> = ({
  mode,
  selectedRegion,
  vcfPath,
  loading,
  publicVcfFolder,
  cohortFolder,
  regions,
  onRegionSelect,
  onOpenSessionManager,
  onClearFilter,
  getSessionData
}) => {
  const [activeTab, setActiveTab] = useState<'visualization' | 'statistics'>('visualization');
  const [showShortcutsHelp, setShowShortcutsHelp] = useState(false);
  const [currentPage, setCurrentPage] = useState<'main' | 'instructions' | 'features'>('main');
  const searchInputRef = useRef<HTMLInputElement>(null);
  const { theme, toggleTheme } = useTheme();

  const shortcuts = createAppShortcuts({
    onSearch: () => {
      searchInputRef.current?.focus();
    },
    onHelp: () => {
      setShowShortcutsHelp(true);
    },
  });

  useKeyboardShortcuts(shortcuts);

  // Render Instructions page
  if (currentPage === 'instructions') {
    return (
      <div className="main-content">
        <div className="top-navbar">
          <div className="navbar-content">
            <div className="navbar-brand">
              <span className="navbar-icon">🧬</span>
              <div className="navbar-title">
                <h1>ProleTRact</h1>
                <span className="navbar-subtitle">TandemTwister visualization tool</span>
              </div>
            </div>
            <div className="navbar-actions">
              <button 
                className="nav-link-btn"
                onClick={() => setCurrentPage('main')}
              >
                ← Back to Analysis
              </button>
              <button 
                className="nav-link-btn"
                onClick={() => setCurrentPage('features')}
              >
                🚀 Features
              </button>
              <button 
                className="theme-toggle-btn"
                onClick={toggleTheme}
                title={theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode'}
              >
                {theme === 'dark' ? '☀️' : '🌙'}
              </button>
            </div>
          </div>
        </div>
        <div className="content-wrapper">
          <InstructionsPage />
        </div>
      </div>
    );
  }

  // Render Features page
  if (currentPage === 'features') {
    return (
      <div className="main-content">
        <div className="top-navbar">
          <div className="navbar-content">
            <div className="navbar-brand">
              <span className="navbar-icon">🧬</span>
              <div className="navbar-title">
                <h1>ProleTRact</h1>
                <span className="navbar-subtitle">TandemTwister visualization tool</span>
              </div>
            </div>
            <div className="navbar-actions">
              <button 
                className="nav-link-btn"
                onClick={() => setCurrentPage('main')}
              >
                ← Back to Analysis
              </button>
              <button 
                className="nav-link-btn"
                onClick={() => setCurrentPage('instructions')}
              >
                📚 Instructions
              </button>
              <button 
                className="theme-toggle-btn"
                onClick={toggleTheme}
                title={theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode'}
              >
                {theme === 'dark' ? '☀️' : '🌙'}
              </button>
            </div>
          </div>
        </div>
        <div className="content-wrapper">
          <FeaturesPage />
        </div>
      </div>
    );
  }

  // Render main analysis page
  return (
    <div className="main-content">
      {/* Top Navigation Bar */}
      <div className="top-navbar">
        <div className="navbar-content">
          <div className="navbar-brand">
            <span className="navbar-icon">🧬</span>
            <div className="navbar-title">
              <h1>ProleTRact</h1>
              <span className="navbar-subtitle">TandemTwister visualization tool</span>
            </div>
          </div>
          <div className="navbar-actions">
            <button 
              className="nav-link-btn"
              onClick={() => setCurrentPage('instructions')}
            >
              📚 Instructions
            </button>
            <button 
              className="nav-link-btn"
              onClick={() => setCurrentPage('features')}
            >
              🚀 Features
            </button>
            <button 
              className="theme-toggle-btn"
              onClick={toggleTheme}
              title={theme === 'dark' ? 'Switch to light mode' : 'Switch to dark mode'}
              aria-label="Toggle theme"
            >
              {theme === 'dark' ? '☀️' : '🌙'}
            </button>
          </div>
        </div>
      </div>
      
      {/* Main Content Area */}
      <div className="content-wrapper">
        {(mode === 'cohort-read' || mode === 'cohort-assembly') ? (
          // Cohort Mode View
          <CohortModeView 
            mode={mode}
            cohortFolder={cohortFolder || ''}
            publicVcfFolder={publicVcfFolder || ''}
            loading={loading}
            selectedRegion={selectedRegion}
            onRegionSelect={onRegionSelect}
          />
        ) : (
          // Individual Mode View
          <>
            {!vcfPath && (
              <div className="welcome-screen">
                <div className="welcome-content">
                  <div className="welcome-logo-video">
                    <video
                      autoPlay
                      loop
                      muted
                      playsInline
                      className="intro-logo-video"
                    >
                      <source src="/intro_logo.mp4" type="video/mp4" />
                    </video>
                    <div className="video-attribution">
                      <a 
                        href="https://iconscout.com/lottie-animations/dna" 
                        target="_blank" 
                        rel="noopener noreferrer"
                        className="attribution-link"
                      >
                        DNA Sequence Cloning And Recombination
                      </a>
                      {' by '}
                      <a 
                        href="https://iconscout.com/contributors/iconsx" 
                        target="_blank" 
                        rel="noopener noreferrer"
                        className="attribution-link"
                      >
                        IconsX
                      </a>
                    </div>
                  </div>
                  <h2>Welcome to ProleTRact</h2>
                  <p>Explore, visualize, and compare tandem repeat regions from TandemTwister outputs</p>
                  
                  {/* Animated Tandem Repeat Explanation */}
                  <TandemRepeatAnimation />
                  
                  <div className="welcome-features">
                    <div className="feature-card">
                      <span className="feature-icon">📁</span>
                      <h3>Load VCF Files</h3>
                      <p>Import your TandemTwister VCF files for analysis</p>
                    </div>
                    <div className="feature-card">
                      <span className="feature-icon">🔬</span>
                      <h3>Region Analysis</h3>
                      <p>Detailed visualization of tandem repeat regions</p>
                    </div>
                    <div className="feature-card">
                      <span className="feature-icon">🌐</span>
                      <h3>Population Comparison</h3>
                      <p>Compare regions across population samples</p>
                    </div>
                    <div className="feature-card">
                      <span className="feature-icon">📊</span>
                      <h3>Statistics Dashboard</h3>
                      <p>Comprehensive statistics and insights</p>
                    </div>
                  </div>
                </div>
              </div>
            )}

        {loading && (
          <div className="loading-container modern">
            <div className="loading-spinner-large">
              <div className="spinner-ring"></div>
              <div className="spinner-ring"></div>
              <div className="spinner-ring"></div>
            </div>
            <span className="loading-text">Loading data...</span>
          </div>
        )}

        {vcfPath && (
          <>
            {/* Modern Tabs */}
            <div className="content-tabs modern">
              <button
                className={`tab-button modern ${activeTab === 'statistics' ? 'active' : ''}`}
                onClick={() => setActiveTab('statistics')}
              >
                <span className="tab-icon">📊</span>
                <span className="tab-label">Statistics</span>
                {activeTab === 'statistics' && <div className="tab-indicator"></div>}
              </button>
              <button
                className={`tab-button modern ${activeTab === 'visualization' ? 'active' : ''}`}
                onClick={() => setActiveTab('visualization')}
              >
                <span className="tab-icon">🔬</span>
                <span className="tab-label">Region Analysis</span>
                {activeTab === 'visualization' && <div className="tab-indicator"></div>}
              </button>
            </div>

            {/* Tab Content - Keep both mounted but hide inactive ones */}
            <div className={`tab-content ${activeTab === 'statistics' ? 'tab-active' : 'tab-hidden'}`}>
              <StatisticsDashboard vcfPath={vcfPath} />
            </div>

            <div className={`tab-content ${activeTab === 'visualization' ? 'tab-active' : 'tab-hidden'}`}>
              {selectedRegion ? (
                <div className="visualization-area">
                  <RegionVisualization 
                    region={selectedRegion} 
                    vcfPath={vcfPath}
                    publicVcfFolder={publicVcfFolder}
                    onRegionSelect={onRegionSelect}
                  />
                </div>
              ) : (
                <div className="no-region-selected">
                  {regions && regions.total_matching === 0 && regions.records.length === 0 && onClearFilter ? (
                    <>
                      <div className="empty-state-icon">🔍</div>
                      <h3>No regions match the current filter</h3>
                      <p>No pathogenic regions were found, or no regions match your filter criteria.</p>
                      <button
                        className="btn btn-primary btn-modern clear-filter-btn"
                        onClick={onClearFilter}
                      >
                        Clear filter &amp; show all regions
                      </button>
                    </>
                  ) : (
                    <>
                      <div className="empty-state-icon">🗺️</div>
                      <h3>No region selected</h3>
                      <p>Please select a region from the navigation below to view its visualization.</p>
                    </>
                  )}
                </div>
              )}
            </div>
          </>
        )}
          </>
        )}
        <BookmarksPanel 
          currentRegion={selectedRegion} 
          onBookmarkSelect={onRegionSelect}
          onOpenSessionManager={onOpenSessionManager}
          getSessionData={getSessionData}
        />
        <KeyboardShortcutsHelp
          shortcuts={shortcuts}
          isOpen={showShortcutsHelp}
          onClose={() => setShowShortcutsHelp(false)}
        />
        <Tutorial />
      </div>
    </div>
  );
};

export default MainContent;

