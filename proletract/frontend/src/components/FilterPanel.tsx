import React from 'react';
import FilterBuilder from './FilterBuilder';
import './FilterPanel.css';

interface FilterPanelProps {
  isOpen: boolean;
  onClose: () => void;
  vcfPath: string;
  availableGenotypes: string[];
  selectedGenotypes: string[];
  availableChromosomes: string[];
  loading: boolean;
  onApplyFilter: (criteria: any) => void;
}

const FilterPanel: React.FC<FilterPanelProps> = ({
  isOpen,
  onClose,
  vcfPath,
  availableGenotypes,
  selectedGenotypes,
  availableChromosomes,
  loading,
  onApplyFilter
}) => {
  if (!isOpen) return null;

  return (
    <>
      <div className="filter-panel-backdrop" onClick={onClose} aria-hidden="true" />
      <div className="filter-panel" role="dialog" aria-label="Advanced filters">
        <div className="filter-panel-header">
          <h2 className="filter-panel-title">
            <span className="filter-panel-icon">🔍</span>
            Filters
          </h2>
          <button
            className="filter-panel-close"
            onClick={onClose}
            aria-label="Close filters"
          >
            ✕
          </button>
        </div>
        <div className="filter-panel-body">
          <FilterBuilder
            vcfPath={vcfPath}
            availableGenotypes={availableGenotypes}
            selectedGenotypes={selectedGenotypes}
            availableChromosomes={availableChromosomes}
            loading={loading}
            onApplyFilter={onApplyFilter}
          />
        </div>
      </div>
    </>
  );
};

export default FilterPanel;
