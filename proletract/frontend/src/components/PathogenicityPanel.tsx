import React, { useState, useEffect } from 'react';
import axios from 'axios';
import './PathogenicityPanel.css';

const API_BASE = process.env.REACT_APP_API_URL || 'http://localhost:8502';

interface PathogenicRegion {
  chr: string;
  start: number;
  end: number;
  gene?: string;
  disease?: string;
  inheritance?: string;
  normal_range?: string;
  pathogenic_threshold?: number;
  description?: string;
  motif?: string;
}

interface PathogenicityPanelProps {
  chr: string;
  pos: number;
  stop: number;
  region: string;
  record?: {
    motif_ids_h1?: string[];
    motif_ids_h2?: string[];
  } | null;
}

const PathogenicityPanel: React.FC<PathogenicityPanelProps> = ({
  chr,
  pos,
  stop,
  region,
  record
}) => {
  const [pathogenicInfo, setPathogenicInfo] = useState<PathogenicRegion | null>(null);
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    // Check if region overlaps with pathogenic catalog
    const checkPathogenicity = async () => {
      setLoading(true);
      setError(null);
      try {
        // This endpoint would need to be implemented in the backend
        const response = await axios.get(`${API_BASE}/api/pathogenic/check`, {
          params: {
            chr: chr,
            start: pos,
            end: stop
          }
        });
        
        // #region agent log
        fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'PathogenicityPanel.tsx:56',message:'Pathogenic check response received',data:{pathogenic:response.data.pathogenic,hasGene:!!response.data.gene,hasDisease:!!response.data.disease,hasInheritance:!!response.data.inheritance,gene:response.data.gene,disease:response.data.disease,inheritance:response.data.inheritance,threshold:response.data.pathogenic_threshold,fullResponse:response.data},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'PATHOGENIC_PANEL'})}).catch(()=>{});
        // #endregion
        
        if (response.data.pathogenic) {
          // #region agent log
          fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'PathogenicityPanel.tsx:62',message:'Setting pathogenicInfo in PathogenicityPanel',data:{gene:response.data.gene,typeOfGene:typeof response.data.gene,disease:response.data.disease,inheritance:response.data.inheritance,threshold:response.data.pathogenic_threshold,fullData:response.data},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'PATHOGENIC_PANEL'})}).catch(()=>{});
          // #endregion
          // Ensure we're setting all fields correctly - preserve values from API response
          // Helper to safely convert and validate string values
          const safeString = (value: any): string | undefined => {
            // Handle null, undefined, or empty values
            if (value === null || value === undefined || value === '') {
              return undefined;
            }
            // Convert to string and trim whitespace
            const str = String(value).trim();
            // Return undefined for empty strings after trimming
            return str.length > 0 ? str : undefined;
          };
          
          // Extract and convert values
          const geneValue = safeString(response.data.gene);
          const diseaseValue = safeString(response.data.disease);
          const inheritanceValue = safeString(response.data.inheritance);
          
          // Debug logging to help identify issues
          console.log('PathogenicityPanel: API Response Data', {
            raw: {
              gene: response.data.gene,
              disease: response.data.disease,
              inheritance: response.data.inheritance,
              threshold: response.data.pathogenic_threshold
            },
            processed: {
              gene: geneValue,
              disease: diseaseValue,
              inheritance: inheritanceValue,
              threshold: response.data.pathogenic_threshold
            },
            types: {
              gene: typeof response.data.gene,
              disease: typeof response.data.disease,
              inheritance: typeof response.data.inheritance
            }
          });
          
          const newPathogenicInfo = {
            chr: response.data.chr,
            start: response.data.start,
            end: response.data.end,
            gene: geneValue,
            disease: diseaseValue,
            inheritance: inheritanceValue,
            pathogenic_threshold: response.data.pathogenic_threshold,
            motif: safeString(response.data.motif),
            normal_range: safeString(response.data.normal_range),
            description: safeString(response.data.description)
          };
          
          console.log('PathogenicityPanel: Setting state with', newPathogenicInfo);
          setPathogenicInfo(newPathogenicInfo);
        } else {
          // #region agent log
          fetch('http://localhost:7245/ingest/c0592157-b6df-40a3-9d0d-4fc90d20aec3',{method:'POST',headers:{'Content-Type':'application/json'},body:JSON.stringify({location:'PathogenicityPanel.tsx:67',message:'Pathogenic check returned false, clearing info',data:{},timestamp:Date.now(),sessionId:'debug-session',runId:'run1',hypothesisId:'PATHOGENIC_PANEL'})}).catch(()=>{});
          // #endregion
          setPathogenicInfo(null);
        }
      } catch (err: any) {
        // Silently fail if endpoint doesn't exist yet
        if (err.response?.status !== 404) {
          console.error('Error checking pathogenicity:', err);
          setError(err.response?.data?.detail || 'Failed to check pathogenicity');
        }
        setPathogenicInfo(null);
      } finally {
        setLoading(false);
      }
    };

    if (chr && pos && stop) {
      checkPathogenicity();
    }
  }, [chr, pos, stop]);

  if (loading) {
    return (
      <div className="pathogenicity-panel loading">
        <div className="pathogenicity-header">
          <h3>🏥 Pathogenicity Information</h3>
        </div>
        <div className="pathogenicity-content">
          <p>Checking pathogenic catalog...</p>
        </div>
      </div>
    );
  }

  if (error) {
    return null; // Don't show panel if there's an error
  }

  if (!pathogenicInfo) {
    return null; // Don't show panel if region is not pathogenic
  }

  // Check if any allele actually exceeds the pathogenic threshold
  const hasPathogenicThreshold = pathogenicInfo.pathogenic_threshold !== undefined && pathogenicInfo.pathogenic_threshold !== null && pathogenicInfo.pathogenic_threshold > 0;
  let hasPathogenicAllele = false;
  
  if (hasPathogenicThreshold && record && pathogenicInfo.pathogenic_threshold !== undefined) {
    // Helper function to count motifs (excluding empty/null/'.' values)
    const calculateMotifCount = (motifIds: string[] | undefined): number => {
      if (!motifIds || !Array.isArray(motifIds)) return 0;
      return motifIds.filter(id => id && id !== '.' && id !== '').length;
    };
    
    const threshold = pathogenicInfo.pathogenic_threshold;
    const h1Count = calculateMotifCount(record.motif_ids_h1);
    const h2Count = calculateMotifCount(record.motif_ids_h2);
    
    // Check if either allele exceeds the threshold
    hasPathogenicAllele = h1Count >= threshold || h2Count >= threshold;
  }

  // Only show as "pathogenic" (red) if there's a threshold AND at least one allele exceeds it
  const isPathogenic = hasPathogenicThreshold && hasPathogenicAllele;
  const copyNumber = Math.floor((stop - pos) / 100); // Rough estimate, would need actual CN from data

  return (
    <div className={`pathogenicity-panel ${isPathogenic ? 'pathogenic' : hasPathogenicThreshold ? 'benign' : 'neutral'}`}>
      <div className="pathogenicity-header">
        <h3>
          {isPathogenic ? '⚠️ Pathogenic Region Detected' : hasPathogenicThreshold ? '✅ Region Below Threshold' : 'ℹ️ Pathogenic Region Information'}
        </h3>
        {isPathogenic && (
          <span className="pathogenicity-badge pathogenic-badge">Pathogenic</span>
        )}
        {hasPathogenicThreshold && !isPathogenic && (
          <span className="pathogenicity-badge benign-badge">Below Threshold</span>
        )}
      </div>
      
      <div className="pathogenicity-content">
        {/* Gene field */}
        {pathogenicInfo.gene && (
          <div className="pathogenicity-item">
            <span className="pathogenicity-label">Gene:</span>
            <span className="pathogenicity-value">{pathogenicInfo.gene}</span>
          </div>
        )}
        
        {/* Disease field */}
        {pathogenicInfo.disease && (
          <div className="pathogenicity-item">
            <span className="pathogenicity-label">Disease:</span>
            <span className="pathogenicity-value">{pathogenicInfo.disease}</span>
          </div>
        )}
        
        {/* Inheritance field */}
        {pathogenicInfo.inheritance && (
          <div className="pathogenicity-item">
            <span className="pathogenicity-label">Inheritance:</span>
            <span className="pathogenicity-value">{pathogenicInfo.inheritance}</span>
          </div>
        )}
        
        {pathogenicInfo.normal_range && (
          <div className="pathogenicity-item">
            <span className="pathogenicity-label">Normal Range:</span>
            <span className="pathogenicity-value">{pathogenicInfo.normal_range}</span>
          </div>
        )}
        
        {pathogenicInfo.pathogenic_threshold !== undefined && (
          <div className="pathogenicity-item">
            <span className="pathogenicity-label">Pathogenic Threshold:</span>
            <span className="pathogenicity-value threshold">
              {pathogenicInfo.pathogenic_threshold} repeats
            </span>
          </div>
        )}
        
        {pathogenicInfo.description && (
          <div className="pathogenicity-description">
            <p>{pathogenicInfo.description}</p>
          </div>
        )}
        
        <div className="pathogenicity-actions">
          <a
            href={`https://www.omim.org/search/?search=${pathogenicInfo.gene || region}`}
            target="_blank"
            rel="noopener noreferrer"
            className="pathogenicity-link"
          >
            View in OMIM →
          </a>
          <a
            href={`https://www.ncbi.nlm.nih.gov/clinvar/?term=${chr}:${pos}`}
            target="_blank"
            rel="noopener noreferrer"
            className="pathogenicity-link"
          >
            View in ClinVar →
          </a>
        </div>
      </div>
    </div>
  );
};

export default PathogenicityPanel;




