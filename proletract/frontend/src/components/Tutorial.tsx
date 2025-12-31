import React, { useState, useEffect } from 'react';
import './Tutorial.css';

interface TutorialStep {
  id: string;
  title: string;
  content: string;
  target?: string; // CSS selector for element to highlight
  position?: 'top' | 'bottom' | 'left' | 'right';
}

const TUTORIAL_STEPS: TutorialStep[] = [
  {
    id: 'welcome',
    title: 'Welcome to ProleTRact!',
    content: 'ProleTRact is a powerful tool for analyzing tandem repeat regions from TandemTwister outputs. This tutorial will guide you through the main features.',
    position: 'bottom',
  },
  {
    id: 'load-vcf',
    title: 'Load VCF File',
    content: 'Start by loading your VCF file. Enter the path to your VCF file in the sidebar and click "Load VCF". The file should be bgzipped and tabix-indexed.',
    target: '.sidebar-section:first-of-type',
    position: 'right',
  },
  {
    id: 'filter-genotypes',
    title: 'Filter by Genotype',
    content: 'Use the genotype filter to focus on specific genotypes. Select or deselect genotypes to filter the regions displayed.',
    target: '.filter-accordion',
    position: 'right',
  },
  {
    id: 'select-region',
    title: 'Select a Region',
    content: 'Choose a region from the navigation panel or use the search bar to find a specific region. Click on a region to view its detailed visualization.',
    target: '.floating-navigation',
    position: 'left',
  },
  {
    id: 'visualization',
    title: 'Region Visualization',
    content: 'The visualization shows the reference and alternate alleles with color-coded motifs. Scroll horizontally to see the full sequence.',
    target: '.region-visualization',
    position: 'top',
  },
  {
    id: 'external-links',
    title: 'External Resources',
    content: 'Quick links to external databases like UCSC Genome Browser, Ensembl, dbSNP, and ClinVar are available for each region.',
    target: '.external-links-panel',
    position: 'top',
  },
  {
    id: 'export',
    title: 'Export Data',
    content: 'Export your visualizations as PNG or SVG, or export data as CSV, JSON, or FASTA format using the export menu.',
    target: '.export-menu-container',
    position: 'bottom',
  },
  {
    id: 'bookmarks',
    title: 'Bookmarks',
    content: 'Bookmark important regions for quick access later. Use the bookmark button in the bottom right corner.',
    target: '.bookmarks-panel',
    position: 'top',
  },
];

interface TutorialProps {
  onComplete?: () => void;
}

const Tutorial: React.FC<TutorialProps> = ({ onComplete }) => {
  const [currentStep, setCurrentStep] = useState(0);
  const [isOpen, setIsOpen] = useState(false);
  const [hasSeenTutorial, setHasSeenTutorial] = useState(false);

  useEffect(() => {
    const seen = localStorage.getItem('proletract_tutorial_seen');
    setHasSeenTutorial(seen === 'true');
    if (!seen) {
      setIsOpen(true);
    }
  }, []);

  const handleNext = () => {
    if (currentStep < TUTORIAL_STEPS.length - 1) {
      setCurrentStep(currentStep + 1);
    } else {
      handleComplete();
    }
  };

  const handlePrevious = () => {
    if (currentStep > 0) {
      setCurrentStep(currentStep - 1);
    }
  };

  const handleSkip = () => {
    handleComplete();
  };

  const handleComplete = () => {
    setIsOpen(false);
    localStorage.setItem('proletract_tutorial_seen', 'true');
    if (onComplete) {
      onComplete();
    }
  };

  const handleStartTutorial = () => {
    setIsOpen(true);
    setCurrentStep(0);
  };

  if (!isOpen) {
    return (
      <button
        className="tutorial-trigger-button"
        onClick={handleStartTutorial}
        title="Start tutorial"
      >
        <span className="tutorial-icon">❓</span>
        <span className="tutorial-label">Tutorial</span>
      </button>
    );
  }

  const step = TUTORIAL_STEPS[currentStep];
  const targetElement = step.target ? document.querySelector(step.target) : null;

  return (
    <>
      {targetElement && (
        <div
          className="tutorial-highlight"
          style={{
            position: 'absolute',
            top: targetElement.getBoundingClientRect().top + window.scrollY,
            left: targetElement.getBoundingClientRect().left + window.scrollX,
            width: targetElement.getBoundingClientRect().width,
            height: targetElement.getBoundingClientRect().height,
          }}
        />
      )}
      <div className="tutorial-overlay" onClick={handleSkip}>
        <div
          className="tutorial-modal"
          onClick={(e) => e.stopPropagation()}
          style={
            targetElement
              ? {
                  position: 'absolute',
                  top:
                    step.position === 'bottom'
                      ? targetElement.getBoundingClientRect().bottom + window.scrollY + 20
                      : step.position === 'top'
                      ? targetElement.getBoundingClientRect().top + window.scrollY - 200
                      : step.position === 'right'
                      ? targetElement.getBoundingClientRect().top + window.scrollY
                      : targetElement.getBoundingClientRect().top + window.scrollY,
                  left:
                    step.position === 'right'
                      ? targetElement.getBoundingClientRect().right + window.scrollX + 20
                      : step.position === 'left'
                      ? targetElement.getBoundingClientRect().left + window.scrollX - 350
                      : targetElement.getBoundingClientRect().left + window.scrollX,
                }
              : {}
          }
        >
          <div className="tutorial-header">
            <div className="tutorial-progress">
              Step {currentStep + 1} of {TUTORIAL_STEPS.length}
            </div>
            <button className="tutorial-close" onClick={handleSkip}>
              ✕
            </button>
          </div>
          <div className="tutorial-content">
            <h3>{step.title}</h3>
            <p>{step.content}</p>
          </div>
          <div className="tutorial-actions">
            <button
              className="tutorial-button tutorial-button-secondary"
              onClick={handlePrevious}
              disabled={currentStep === 0}
            >
              Previous
            </button>
            <button
              className="tutorial-button tutorial-button-skip"
              onClick={handleSkip}
            >
              Skip Tutorial
            </button>
            <button
              className="tutorial-button tutorial-button-primary"
              onClick={handleNext}
            >
              {currentStep === TUTORIAL_STEPS.length - 1 ? 'Finish' : 'Next'}
            </button>
          </div>
        </div>
      </div>
    </>
  );
};

export default Tutorial;




