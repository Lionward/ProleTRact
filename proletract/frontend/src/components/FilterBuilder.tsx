import React, { useState, useEffect } from 'react';
import './FilterBuilder.css';

interface FilterCriteria {
  motifSizeMin?: number;
  motifSizeMax?: number;
  cnMin?: number;
  cnMax?: number;
  chromosomes?: string[];
  genotypes?: string[];
  pathogenicOnly?: boolean;
  hasAnnotations?: boolean;
  annotationTags?: string[];
}

interface FilterPreset {
  id: string;
  name: string;
  criteria: FilterCriteria;
  createdAt: number;
}

interface FilterBuilderProps {
  vcfPath: string;
  availableGenotypes: string[];
  selectedGenotypes?: string[];
  availableChromosomes?: string[];
  loading?: boolean;
  onFilterChange?: (criteria: FilterCriteria) => void;
  onApplyFilter?: (criteria: FilterCriteria) => void;
}

interface QuickFilterDef {
  name: string;
  criteria: FilterCriteria;
  isActive: (c: FilterCriteria) => boolean;
  clearCriteria: Partial<FilterCriteria>;
}

const QUICK_FILTERS: QuickFilterDef[] = [
  { name: 'High CN (≥10)', criteria: { cnMin: 10 }, isActive: (c) => c.cnMin === 10, clearCriteria: { cnMin: undefined } },
  { name: 'Pathogenic', criteria: { pathogenicOnly: true }, isActive: (c) => c.pathogenicOnly === true, clearCriteria: { pathogenicOnly: false } }
];

const FilterBuilder: React.FC<FilterBuilderProps> = ({
  vcfPath,
  availableGenotypes,
  selectedGenotypes = availableGenotypes,
  availableChromosomes = [],
  loading = false,
  onFilterChange,
  onApplyFilter
}) => {
  // Load criteria from localStorage for persistence across panel open/close
  const [criteria, setCriteria] = useState<FilterCriteria>(() => {
    const stored = localStorage.getItem('filter_criteria_current');
    if (stored) {
      try {
        return JSON.parse(stored) as FilterCriteria;
      } catch (err) {
        return {};
      }
    }
    return {};
  });
  const [presets, setPresets] = useState<FilterPreset[]>([]);
  const [presetName, setPresetName] = useState('');
  const [showPresetDialog, setShowPresetDialog] = useState(false);
  const [selectedPreset, setSelectedPreset] = useState<string | null>(null);

  // Load presets from localStorage
  useEffect(() => {
    const stored = localStorage.getItem('filter_presets');
    if (stored) {
      try {
        const parsed = JSON.parse(stored) as FilterPreset[];
        setPresets(parsed);
      } catch (err) {
        console.error('Failed to load filter presets:', err);
      }
    }
  }, []);

  // Save criteria to localStorage whenever it changes
  useEffect(() => {
    localStorage.setItem('filter_criteria_current', JSON.stringify(criteria));
  }, [criteria]);

  // Save presets to localStorage
  const savePresets = (newPresets: FilterPreset[]) => {
    localStorage.setItem('filter_presets', JSON.stringify(newPresets));
    setPresets(newPresets);
  };

  // Update criteria
  const updateCriteria = (updates: Partial<FilterCriteria>) => {
    const newCriteria = { ...criteria, ...updates };
    setCriteria(newCriteria);
    if (onFilterChange) {
      onFilterChange(newCriteria);
    }
  };

  // Toggle quick filter (press/unpress like genotype buttons)
  const toggleQuickFilter = (quick: QuickFilterDef) => {
    const isActive = quick.isActive(criteria);
    const newCriteria = isActive
      ? { ...criteria, ...quick.clearCriteria }
      : { ...criteria, ...quick.criteria };
    setCriteria(newCriteria);
    if (onApplyFilter) {
      onApplyFilter(newCriteria);
    }
  };

  // Save preset
  const savePreset = () => {
    if (!presetName.trim()) {
      alert('Please enter a name for the preset');
      return;
    }

    const newPreset: FilterPreset = {
      id: `preset_${Date.now()}`,
      name: presetName.trim(),
      criteria: { ...criteria },
      createdAt: Date.now()
    };

    const newPresets = [...presets, newPreset];
    savePresets(newPresets);
    setPresetName('');
    setShowPresetDialog(false);
  };

  // Load preset
  const loadPreset = (presetId: string) => {
    const preset = presets.find(p => p.id === presetId);
    if (preset) {
      setCriteria(preset.criteria);
      setSelectedPreset(presetId);
      if (onApplyFilter) {
        onApplyFilter(preset.criteria);
      }
    }
  };

  // Delete preset
  const deletePreset = (presetId: string) => {
    if (window.confirm('Are you sure you want to delete this preset?')) {
      const newPresets = presets.filter(p => p.id !== presetId);
      savePresets(newPresets);
      if (selectedPreset === presetId) {
        setSelectedPreset(null);
      }
    }
  };

  // Clear filter
  const clearFilter = () => {
    setCriteria({});
    localStorage.removeItem('filter_criteria_current');
    setSelectedPreset(null);
    if (onApplyFilter) {
      onApplyFilter({});
    }
  };

  // Toggle chromosome
  const toggleChromosome = (chr: string) => {
    const current = criteria.chromosomes || [];
    const updated = current.includes(chr)
      ? current.filter(c => c !== chr)
      : [...current, chr];
    updateCriteria({ chromosomes: updated.length > 0 ? updated : undefined });
  };

  // Toggle genotype
  const toggleGenotype = (gt: string) => {
    const current = criteria.genotypes || [];
    const updated = current.includes(gt)
      ? current.filter(g => g !== gt)
      : [...current, gt];
    updateCriteria({ genotypes: updated.length > 0 ? updated : undefined });
  };

  // Apply filter - use selectedGenotypes when criteria.genotypes is undefined
  const applyFilter = () => {
    const payload = { ...criteria, genotypes: criteria.genotypes ?? selectedGenotypes };
    if (onApplyFilter) {
      onApplyFilter(payload);
    }
  };

  return (
    <div className="filter-builder">
      <div className="filter-builder-header">
        <button
          className="filter-btn filter-btn-secondary"
          onClick={clearFilter}
          title="Clear all filters"
        >
          Clear
        </button>
      </div>

      {/* Quick Filters - toggle buttons like genotype chips */}
      <div className="filter-section">
        <label className="filter-label">Quick Filters</label>
        <div className="quick-filters">
          {QUICK_FILTERS.map((quick, idx) => (
            <button
              key={idx}
              className={`quick-filter-btn ${quick.isActive(criteria) ? 'selected' : ''}`}
              onClick={() => toggleQuickFilter(quick)}
              title={quick.name}
            >
              {quick.name}
            </button>
          ))}
        </div>
      </div>

      {/* Motif Size Range */}
      <div className="filter-section">
        <label className="filter-label">Motif Size Range</label>
        <div className="filter-range">
          <input
            type="number"
            className="filter-input filter-input-number"
            placeholder="Min"
            value={criteria.motifSizeMin || ''}
            onChange={(e) =>
              updateCriteria({
                motifSizeMin: e.target.value ? parseInt(e.target.value, 10) : undefined
              })
            }
          />
          <span className="filter-range-separator">to</span>
          <input
            type="number"
            className="filter-input filter-input-number"
            placeholder="Max"
            value={criteria.motifSizeMax || ''}
            onChange={(e) =>
              updateCriteria({
                motifSizeMax: e.target.value ? parseInt(e.target.value, 10) : undefined
              })
            }
          />
        </div>
      </div>

      {/* Copy Number Range */}
      <div className="filter-section">
        <label className="filter-label">Copy Number Range</label>
        <div className="filter-range">
          <input
            type="number"
            className="filter-input filter-input-number"
            placeholder="Min CN"
            value={criteria.cnMin || ''}
            onChange={(e) =>
              updateCriteria({
                cnMin: e.target.value ? parseFloat(e.target.value) : undefined
              })
            }
          />
          <span className="filter-range-separator">to</span>
          <input
            type="number"
            className="filter-input filter-input-number"
            placeholder="Max CN"
            value={criteria.cnMax || ''}
            onChange={(e) =>
              updateCriteria({
                cnMax: e.target.value ? parseFloat(e.target.value) : undefined
              })
            }
          />
        </div>
      </div>

      {/* Chromosomes */}
      {availableChromosomes.length > 0 && (
        <div className="filter-section">
          <label className="filter-label">Chromosomes</label>
          <div className="filter-chips">
            {availableChromosomes.map(chr => (
              <button
                key={chr}
                className={`filter-chip ${(criteria.chromosomes || []).includes(chr) ? 'selected' : ''}`}
                onClick={() => toggleChromosome(chr)}
              >
                {chr}
              </button>
            ))}
          </div>
        </div>
      )}

      {/* Genotypes - sync with selectedGenotypes when criteria.genotypes is undefined */}
      {availableGenotypes.length > 0 && (
        <div className="filter-section">
          <label className="filter-label">Genotypes</label>
          <div className="filter-chips">
            {availableGenotypes.map(gt => {
              const activeGenotypes = criteria.genotypes ?? selectedGenotypes;
              const isSelected = activeGenotypes.includes(gt);
              return (
                <button
                  key={gt}
                  className={`filter-chip ${isSelected ? 'selected' : ''}`}
                  onClick={() => toggleGenotype(gt)}
                >
                  {gt}
                </button>
              );
            })}
          </div>
        </div>
      )}

      {/* Annotations filter */}
      <div className="filter-section">
        <label className="filter-label">Additional</label>
        <div className="filter-checkboxes">
          <label className="filter-checkbox">
            <input
              type="checkbox"
              checked={criteria.hasAnnotations || false}
              onChange={(e) => updateCriteria({ hasAnnotations: e.target.checked || undefined })}
            />
            <span>Has annotations</span>
          </label>
        </div>
      </div>

      {/* Presets */}
      <div className="filter-section">
        <div className="filter-presets-header">
          <label className="filter-label">Saved Presets</label>
          <button
            className="filter-btn filter-btn-link"
            onClick={() => setShowPresetDialog(true)}
            title="Save current filter as preset"
          >
            + Save Preset
          </button>
        </div>
        {presets.length > 0 ? (
          <div className="filter-presets">
            {presets.map(preset => (
              <div key={preset.id} className="filter-preset-item">
                <button
                  className={`filter-preset-btn ${selectedPreset === preset.id ? 'active' : ''}`}
                  onClick={() => loadPreset(preset.id)}
                  title={preset.name}
                >
                  {preset.name}
                </button>
                <button
                  className="filter-preset-delete"
                  onClick={() => deletePreset(preset.id)}
                  title="Delete preset"
                >
                  ×
                </button>
              </div>
            ))}
          </div>
        ) : (
          <div className="filter-presets-empty">No saved presets</div>
        )}
      </div>

      {/* Apply Button */}
      <div className="filter-actions">
        <button
          className="filter-btn filter-btn-primary"
          onClick={applyFilter}
          disabled={loading}
        >
          {loading ? 'Applying...' : 'Apply Filter'}
        </button>
      </div>

      {/* Save Preset Dialog */}
      {showPresetDialog && (
        <div className="filter-dialog-overlay" onClick={() => setShowPresetDialog(false)}>
          <div className="filter-dialog" onClick={(e) => e.stopPropagation()}>
            <h4 className="filter-dialog-title">Save Filter Preset</h4>
            <input
              type="text"
              className="filter-input"
              placeholder="Preset name"
              value={presetName}
              onChange={(e) => setPresetName(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === 'Enter') {
                  savePreset();
                } else if (e.key === 'Escape') {
                  setShowPresetDialog(false);
                }
              }}
              autoFocus
            />
            <div className="filter-dialog-actions">
              <button
                className="filter-btn filter-btn-secondary"
                onClick={() => setShowPresetDialog(false)}
              >
                Cancel
              </button>
              <button
                className="filter-btn filter-btn-primary"
                onClick={savePreset}
              >
                Save
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default FilterBuilder;
