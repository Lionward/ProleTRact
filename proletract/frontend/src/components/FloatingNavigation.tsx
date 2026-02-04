import React, { useRef, useState, useEffect } from 'react';
import { FilterResponse } from '../types';
import './FloatingNavigation.css';

interface FloatingNavigationProps {
  regions: FilterResponse | null;
  currentPage: number;
  onPageChange: (page: number) => void;
  selectedRegion: string;
  onNextRegion: () => void;
  onPreviousRegion: () => void;
  onJumpToRegion: (regionNumber: number) => void;
  loading: boolean;
}

const FloatingNavigation: React.FC<FloatingNavigationProps> = ({
  regions,
  currentPage,
  onPageChange,
  selectedRegion,
  onNextRegion,
  onPreviousRegion,
  onJumpToRegion,
  loading
}) => {
  const [editingRegion, setEditingRegion] = useState(false);
  const [regionInput, setRegionInput] = useState('');
  const inputRef = useRef<HTMLInputElement>(null);

  // Handle input focus - must be before early return
  useEffect(() => {
    if (editingRegion && inputRef.current) {
      inputRef.current.focus();
      inputRef.current.select();
    }
  }, [editingRegion]);

  // Show navigation if there are regions
  // Hide only if regions is null, or if there's exactly 1 page AND only 1 record (to avoid clutter)
  if (!regions) {
    return null;
  }
  // Only hide if there's exactly 1 page, 1 record, AND total matching is 1 or less
  if (regions.total_pages <= 1 && regions.records.length <= 1 && regions.total_matching <= 1) {
    return null;
  }

  // Check if we can navigate to next/previous region
  const currentIndex = regions.records.findIndex(r => r.region === selectedRegion);
  const canGoNextRegion = currentIndex < regions.records.length - 1 || currentPage < regions.total_pages - 1;
  const canGoPreviousRegion = currentIndex > 0 || currentPage > 0;
  const isFirstRegion = currentIndex === 0 && currentPage === 0;
  const isLastRegion = currentIndex === regions.records.length - 1 && currentPage === regions.total_pages - 1;

  // Calculate current region number (1-indexed)
  const currentRegionNumber = currentIndex >= 0 
    ? (currentPage * 50) + currentIndex + 1 
    : (currentPage * 50) + 1;
  const totalRegions = regions.total_matching;

  // Handle input submit
  const handleRegionInputSubmit = () => {
    const regionNum = parseInt(regionInput.replace(/,/g, ''), 10);
    if (!isNaN(regionNum) && regionNum >= 1 && regionNum <= totalRegions) {
      onJumpToRegion(regionNum);
      setEditingRegion(false);
      setRegionInput('');
    } else {
      alert(`Please enter a number between 1 and ${totalRegions.toLocaleString()}`);
    }
  };

  // Handle input key press
  const handleRegionInputKeyDown = (e: React.KeyboardEvent<HTMLInputElement>) => {
    if (e.key === 'Enter') {
      handleRegionInputSubmit();
    } else if (e.key === 'Escape') {
      setEditingRegion(false);
      setRegionInput('');
    }
  };

  return (
    <div className="floating-navigation">
      {/* Page navigation buttons (on the sides) - using original prev/next classes */}
      <button
        className="floating-nav-btn floating-nav-prev"
        onClick={() => onPageChange(currentPage - 1)}
        disabled={currentPage === 0 || loading}
        title="Previous page (50 regions)"
        aria-label="Previous page"
      >
        <span className="floating-nav-icon">|&lt;&lt;</span>
      </button>
      
      <div className="floating-nav-divider"></div>
      
      {/* Single region navigation buttons (in the middle) */}
      <button
        className="floating-nav-btn floating-nav-prev-region"
        onClick={onPreviousRegion}
        disabled={isFirstRegion || loading}
        title="Previous region (single step)"
        aria-label="Previous region"
      >
        <span className="floating-nav-icon">◀</span>
      </button>
      
      <div className="floating-nav-info">
        {editingRegion ? (
          <input
            ref={inputRef}
            type="text"
            className="floating-nav-input"
            value={regionInput}
            onChange={(e) => setRegionInput(e.target.value)}
            onBlur={handleRegionInputSubmit}
            onKeyDown={handleRegionInputKeyDown}
            placeholder={currentRegionNumber.toString()}
          />
        ) : (
          <>
            <span 
              className="floating-nav-current"
              onClick={() => {
                setEditingRegion(true);
                setRegionInput(currentRegionNumber.toString());
              }}
              title="Click to jump to a specific region number"
              style={{ cursor: 'pointer' }}
            >
              {currentRegionNumber.toLocaleString()}
            </span>
            <span className="floating-nav-separator">/</span>
            <span className="floating-nav-total">{totalRegions.toLocaleString()}</span>
          </>
        )}
      </div>
      
      <button
        className="floating-nav-btn floating-nav-next-region"
        onClick={onNextRegion}
        disabled={isLastRegion || loading}
        title="Next region (single step)"
        aria-label="Next region"
      >
        <span className="floating-nav-icon">▶</span>
      </button>
      
      <div className="floating-nav-divider"></div>
      
      {/* Page navigation buttons (on the sides) - using original prev/next classes */}
      <button
        className="floating-nav-btn floating-nav-next"
        onClick={() => onPageChange(currentPage + 1)}
        disabled={currentPage >= regions.total_pages - 1 || loading}
        title="Next page (50 regions)"
        aria-label="Next page"
      >
        <span className="floating-nav-icon">&gt;&gt;|</span>
      </button>
    </div>
  );
};

export default FloatingNavigation;

