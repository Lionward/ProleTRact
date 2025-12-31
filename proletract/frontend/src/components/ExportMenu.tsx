import React, { useState, useRef } from 'react';
import './ExportMenu.css';
import { exportElementAsPNG, exportElementAsSVG, generateFilename } from '../utils/exportUtils';

interface ExportMenuProps {
  elementRef?: React.RefObject<HTMLElement>;
  elementId?: string;
  filename?: string;
  onExportCSV?: () => void;
  onExportJSON?: () => void;
  onExportFASTA?: () => void;
  showImageExport?: boolean;
  showDataExport?: boolean;
}

const ExportMenu: React.FC<ExportMenuProps> = ({
  elementRef,
  elementId,
  filename = 'export',
  onExportCSV,
  onExportJSON,
  onExportFASTA,
  showImageExport = true,
  showDataExport = true,
}) => {
  const [isOpen, setIsOpen] = useState(false);
  const menuRef = useRef<HTMLDivElement>(null);

  const handleExportPNG = async () => {
    try {
      let element: HTMLElement | null = null;
      
      if (elementRef?.current) {
        element = elementRef.current;
      } else if (elementId) {
        element = document.getElementById(elementId);
      }
      
      if (!element) {
        console.error('Export: Element not found', { elementRef: elementRef?.current, elementId });
        alert('Element not found for export. Please ensure the component is fully loaded.');
        return;
      }
      
      console.log('Exporting PNG:', { element, filename: generateFilename(filename, 'png') });
      await exportElementAsPNG(element, generateFilename(filename, 'png'));
      setIsOpen(false);
    } catch (error) {
      console.error('Error in handleExportPNG:', error);
      alert('Failed to export PNG. Check console for details.');
    }
  };

  const handleExportSVG = () => {
    try {
      let element: HTMLElement | null = null;
      
      if (elementRef?.current) {
        element = elementRef.current;
      } else if (elementId) {
        element = document.getElementById(elementId);
      }
      
      if (!element) {
        console.error('Export: Element not found', { elementRef: elementRef?.current, elementId });
        alert('Element not found for export. Please ensure the component is fully loaded.');
        return;
      }
      
      console.log('Exporting SVG:', { element, filename: generateFilename(filename, 'svg') });
      exportElementAsSVG(element, generateFilename(filename, 'svg'));
      setIsOpen(false);
    } catch (error) {
      console.error('Error in handleExportSVG:', error);
      alert('Failed to export SVG. Check console for details.');
    }
  };

  // Close menu when clicking outside
  React.useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (menuRef.current && !menuRef.current.contains(event.target as Node)) {
        setIsOpen(false);
      }
    };

    if (isOpen) {
      document.addEventListener('mousedown', handleClickOutside);
      return () => {
        document.removeEventListener('mousedown', handleClickOutside);
      };
    }
  }, [isOpen]);

  return (
    <div className="export-menu-container" ref={menuRef}>
      <button
        className="export-menu-button"
        onClick={() => setIsOpen(!isOpen)}
        title="Export options"
      >
        <span className="export-icon">📥</span>
        <span className="export-label">Export</span>
        <span className={`export-arrow ${isOpen ? 'open' : ''}`}>▼</span>
      </button>
      
      {isOpen && (
        <div className="export-menu-dropdown">
          {showImageExport && (
            <>
              <button
                className="export-menu-item"
                onClick={handleExportPNG}
              >
                <span className="export-item-icon">🖼️</span>
                <span>Export as PNG</span>
              </button>
              <button
                className="export-menu-item"
                onClick={handleExportSVG}
              >
                <span className="export-item-icon">📐</span>
                <span>Export as SVG</span>
              </button>
              {showDataExport && <div className="export-menu-divider" />}
            </>
          )}
          
          {showDataExport && (
            <>
              {onExportCSV && (
                <button
                  className="export-menu-item"
                  onClick={() => {
                    try {
                      onExportCSV();
                      setIsOpen(false);
                    } catch (error) {
                      console.error('Error exporting CSV:', error);
                      alert('Failed to export CSV. Check console for details.');
                    }
                  }}
                >
                  <span className="export-item-icon">📊</span>
                  <span>Export as CSV</span>
                </button>
              )}
              {onExportJSON && (
                <button
                  className="export-menu-item"
                  onClick={() => {
                    try {
                      onExportJSON();
                      setIsOpen(false);
                    } catch (error) {
                      console.error('Error exporting JSON:', error);
                      alert('Failed to export JSON. Check console for details.');
                    }
                  }}
                >
                  <span className="export-item-icon">📄</span>
                  <span>Export as JSON</span>
                </button>
              )}
              {onExportFASTA && (
                <button
                  className="export-menu-item"
                  onClick={() => {
                    try {
                      onExportFASTA();
                      setIsOpen(false);
                    } catch (error) {
                      console.error('Error exporting FASTA:', error);
                      alert('Failed to export FASTA. Check console for details.');
                    }
                  }}
                >
                  <span className="export-item-icon">🧬</span>
                  <span>Export as FASTA</span>
                </button>
              )}
            </>
          )}
        </div>
      )}
    </div>
  );
};

export default ExportMenu;


